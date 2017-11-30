function [epis,finalepisize,validvol,meanvol,additionalvol] = preprocessfmri(figuredir,inplanes,inplanesizes, ...
  fieldmaps,fieldmapbrains,fieldmapsizes,fieldmapdeltate,fieldmapunwrap,fieldmapsmoothing, ...
  epis,episize,epiinplanematrixsize,epitr,episliceorder,epiphasedir,epireadouttime,epifieldmapasst, ...
  numepiignore,motionreference,motioncutoff,extratrans,targetres,sliceshiftband, ...
  fmriqualityparams,fieldmaptimeinterp,mcmask,maskoutnans,epiignoremcvol,dformat,epismoothfwhm,wantpushalt)

% function [epis,finalepisize,validvol,meanvol,additionalvol = preprocessfmri(figuredir,inplanes,inplanesizes, ...
%   fieldmaps,fieldmapbrains,fieldmapsizes,fieldmapdeltate,fieldmapsmoothing, ...
%   epis,episize,epiinplanematrixsize,epitr,episliceorder,epiphasedir,epireadouttime,epifieldmapasst, ...
%   numepiignore,motionreference,motioncutoff,extratrans,targetres,sliceshiftband, ...
%   fmriqualityparams,fieldmaptimeinterp,mcmask,maskoutnans,epiignoremcvol,dformat,epismoothfwhm,wantpushalt)
%
% <figuredir> is the directory to write figures and results to.
%   [] means do not write figures and results.
% <inplanes> is a 3D volume or a cell vector of 3D volumes.  can be {}.
% <inplanesizes> is a 3-element vector with the voxel size in mm
%   or a cell vector of such vectors.  should mirror <inplanes>.
%   if a single vector, we automatically repeat that vector for multiple
%   <inplanes>.
% <fieldmaps> is a 3D volume of phase values (in [-pi,pi]) or a cell vector
%   of such volumes or {A B} where A is a cell vector of one or more such
%   volumes and B is a vector of time values to associate with the volumes. 
%   if B is omitted or given as [], we default to 1:N where N is the number of
%   volumes supplied.  <fieldmaps> can be {}, which indicates that no 
%   undistortion is to be applied.
% <fieldmapbrains> is a 3D volume of magnitude brains or a cell vector
%   of such volumes.  can be {}.  should mirror <fieldmaps>.
% <fieldmapsizes> is a 3-element vector with the voxel size in mm
%   or a cell vector of such vectors.  can be {}.  should mirror <fieldmaps>.
%   if a single vector, we automatically repeat that vector for multiple
%   <fieldmaps>.
% <fieldmapdeltate> is the difference in TE that was used for the fieldmap
%   or a vector of such differences.  should be in milliseconds.  can be [].
%   should mirror <fieldmaps>.  if a single value, we automatically repeat that 
%   value for multiple <fieldmaps>.
% <fieldmapunwrap> is whether to attempt to unwrap the fieldmap.  can be:
%   (1) 0 means no.
%   (2) 1 means yes and use the '-s -t 0' flags in prelude.
%   (3) a string with the specific flags to use (e.g. '', '-f').
%   can be a cell vector of things like (1)-(3), and can be {}.  should mirror <fieldmaps>.  
%   if a single element, we automatically repeat that element for multiple <fieldmaps>.
% <fieldmapsmoothing> is a 3-element vector with the size of the
%   bandwidth in mm to use in the local linear regression or a cell vector
%   of such vectors.  can be {}.  should mirror <fieldmaps>.  if a single vector, we 
%   automatically repeat that vector for multiple <fieldmaps>.  the bandwidth 
%   along each dimension must be greater than the voxel size of the fieldmap 
%   along that dimension (this is because we need more than one data point along each 
%   dimension in order to perform local linear regression).  also, bandwidths
%   can be Inf (this results in pure linear regression along the corresponding 
%   dimension).  also, instead of a 3-element vector, you can supply NaN 
%   which means to omit smoothing.  for details on the meaning of the bandwidth,
%   see localregression3d.m, but here is an example: if <fieldmapsmoothing> is
%   [7.5 7.5 7.5], this means that the smoothing kernel should fall to 0 when
%   you are 7.5 mm away from the center point.
% <epis> is a 4D volume or a cell vector of 4D volumes.
%   these volumes should be double format but should be suitable for 
%   interpretation as int16.  there must be at least one EPI run.  
%   the first three dimensions must be consistent across cases.
%   a special case is that <epis> can be phase angles. to specify this
%   case, <epis> must be converted from angles to unit-length complex numbers,
%   multiplied by 10,000, and cast to int16; and you must also supply
%   <wantpushalt>.  when we are all done processing, the output in this special 
%   case is returned as angles expressed as int16 in the range [0,4095], with the
%   intention of being interpreted as [0,2*pi]. some special requirements of the phase
%   angle case: if <episliceorder> is used, it must be of the case where it is
%   a cell vector with at least two elements. also, in the phase data case,
%   we make no guarantees on the validity of anything related to temporalsnr, 
%   <meanvol> and <additionalvol> (so distrust those quantities!).  in fact,
%   <meanvol> is returned as [] and <additionalvol> is returned as {[] []}.
% <episize> is a 3-element vector with the voxel size in mm.
% <epiinplanematrixsize> is [A B] where A and B are the in-plane frequency-encode
%   and phase-encode matrix sizes, respectively.  for example, [76 64] indicates
%   76 frequency-encode steps and 64 phase-encode steps.  the reason that this
%   input is necessary is that the dimensions of <epis> may be larger than
%   the actual measured matrix size due to zero-padding in the reconstruction process.
%   can be [] in which case we default to the size of the first two dimensions 
%   of the EPI data.
% <epitr> is the TR in seconds or a vector of TRs.  should mirror <epis>.
%   if a single TR, we automatically repeat that integer for multiple <epis>.
% <episliceorder> is a vector of positive integers (some permutation
%   of 1:S where S is the number of slices in the EPI volumes).  you can 
%   also specify 'sequential' or 'interleaved' or 'interleavedalt' ---
%   'sequential' translates into 1:S, 'interleaved' translates into 
%   [1:2:S 2:2:S], and 'interleavedalt' translates into [2:2:S 1:2:S]
%   when S is even and [1:2:S 2:2:S] when S is odd.  can also be {X} where 
%   X is a vector of positive integers that specify the times at which each 
%   slice is collected.  for example, if X is [1 2 3 1 2 3] this means that 
%   the 1st and 4th slices were acquired first, then the 2nd and 5th slices, 
%   and then the 3rd and 6th slices.  can also be {X NEWTR} which is the same
%   as the {X} case except that we prepare the data at a new TR (NEWTR) and
%   in doing so we use pchip interpolation (and first and last data point padding)
%   (instead of the usual sinc interpolation).  NEWTR can be a vector in which
%   case different runs get different TRs.  can also be {X NEWTR NEWOFFSET} which
%   allows temporal offsets (positive means start in the future) in seconds, where
%   NEWOFFSET can be a scalar or a vector (specifying different values for different
%   runs).  you can also set <episliceorder> to [] which means do not perform slice 
%   time correction.
% <epiphasedir> is an integer indicating the phase-encode direction or
%   a vector of such integers.  should mirror <epis>.  if a single integer,
%   we automatically repeat that integer for multiple <epis>.
%   this input matters only if <fieldmaps> and/or <sliceshiftband> are provided.
%   valid values are 1, -1, 2, or -2.  a value of 1 means the phase-encode 
%   direction is oriented along the first matrix dimension (e.g. 1->64); 
%   -1 means the reverse (e.g. 64->1).  a value of 2 means the phase-encode
%   direction is oriented along the second matrix dimension (e.g. 1->64);
%   -2 means the reverse (e.g. 64->1).
% <epireadouttime> is the duration of the EPI readout in milliseconds.
%   matters only if fieldmap(s) are provided.
% <epifieldmapasst> indicates which fieldmaps to apply to which EPI runs.  
%   matters only if fieldmap(s) are provided.  if NaN, that means do not 
%   apply any fieldmaps to any run.  should be a cell vector of the 
%   same length as <epis>.  each element should either be a non-negative
%   integer or a sorted two-element vector [G H].  for a given element, if 
%   it is a positive integer, that indicates the index of the fieldmap to use; 
%   if it is 0, that indicates do not apply a fieldmap for that EPI run; if
%   it is a two-element vector, that indicates to interpolate from time value
%   G to time value H in order to figure out the fieldmap values to use
%   for that EPI run (in this case, each volume within the run gets a different
%   fieldmap correction).  if you pass in a vector of non-negative integers,
%   we automatically convert that to a cell vector.  there are some special
%   cases for your convenience.  if <epifieldmapasst> is passed in as [], 
%   that means (1) if there is one fieldmap then apply that fieldmap to all 
%   EPI runs, (2) if there are as many fieldmaps as EPI runs then apply fieldmaps
%   1-to-1 to the EPI runs, (3) if there is one more fieldmap than EPI runs, then
%   interpolate between each successive fieldmap to obtain the fieldmap to use
%   for each EPI run, and (4) otherwise, give an error.
% <numepiignore> (optional) is a non-negative integer indicating number of volumes
%   at the beginning of the EPI run to ignore or a vector of such integers.  
%   should mirror <epis>.  if a single integer, we automatically repeat that 
%   integer for multiple <epis>.  default: 0.  can also be {[A1 A2] [B1 B2] ... [N1 N2]}
%   where the length of this cell vector is the same as the number of <epis>,
%   and where A1, B1, ..., N1 are non-negative integers indicating number of
%   volumes at the beginning of the EPI runs to ignore, and A2, B2, ..., N2 are
%   non-negative integers indicating number of volumes at the end of the EPI runs
%   to ignore.
% <motionreference> (optional) indicates which EPI volume should be used 
%   as the reference for motion correction.  the format is [g h] where g indicates 
%   the index of the EPI run and h indicates the index of the volume within 
%   that run (after ignoring volumes according to <numepiignore>).  can also be the 
%   3D volume itself (we use this case if ~isvector(<motionreference>)).
%   can also be NaN which indicates that no motion correction should be performed.
%   default: [1 1].
% <motioncutoff> (optional) is the low-pass filter cutoff in Hz to use for
%   motion parameter estimates.  can be Inf which in effect means to not
%   perform any low-pass filtering.  this input matters only if motion 
%   correction is performed.  default: 1/90.
% <extratrans> (optional) is a 4x4 transformation matrix that maps points in the
%   matrix space of the EPI volumes to a new location.  if supplied, then volumes 
%   will be resampled at the new location.  for example, if <extratrans>
%   is [1 0 0 1; 0 1 0 0; 0 0 1 0; 0 0 0 1], then this will cause volumes to be 
%   resampled at a location corresponding to a one-voxel shift along the first
%   dimension.  <extratrans> can also be {X} where X is 4 x vertices, indicating
%   the exact locations (relative to the matrix space of the 3D volume) at which
%   to sample the data.  in this case, <targetres> must be [] and <finalepisize>
%   is returned as [].  default: eye(4).
% <targetres> (optional) is
%   (1) [X Y Z], a 3-element vector with the number of voxels desired for the final 
%       output.  if supplied, then volumes will be interpolated only at the points 
%       necessary to achieve the desired <targetres>.  we assume that the field-
%       of-view is to be preserved.
%   (2) {[A B C] [D E F] G H}, where [A B C] is a 3-element vector with the number of voxels
%       desired for the final output and [D E F] is the voxel size in mm.  in this case, 
%       we do not assume the field-of-view is preserved; rather, we just assume the 
%       desired grid is ndgrid(1:A,1:B,1:C) (after potentially moving according to
%       <extratrans>).  G is 0/1 indicating whether to tightly crop the output 
%       volumes (see below for more details).  H is 0/1 indicating whether to 
%       save only the <validvol> voxels, and we do so in a flattened format 
%       (A x 1 x 1 x T).
%   default is [] which means to do nothing special.
% <sliceshiftband> (optional) is [LOW HIGH] where LOW and HIGH are frequencies in Hz
%   defining the range of a band-pass filter.  if supplied, we perform the slice-shifting
%   correction based on band-pass filtering the center-of-mass calculations according
%   to [LOW HIGH].  default is [] which means to omit the slice-shifting correction.
%   reasonable values for LOW and HIGH might be 1/360 and 1/20, respectively.
%   LOW can be 0 which means only low-pass filter with a cutoff of HIGH.
% <fmriqualityparams> (optional) is {seed numrandom knobs}, as described in
%   fmriquality.m.  default is [] which means to use the defaults in fmriquality.m.
%   special case is NaN which means to skip the fmriquality stuff.
% <fieldmaptimeinterp> (optional) is 'linear' or 'cubic' indicating the type of
%   interpolation to use for the [G H] cases in <epifieldmapasst>.  note that
%   to perform interpolation, we use the interp1.m function and in particular we 
%   use 'extrap' (thus, it is possible to use values for G and H that are outside 
%   the range of the original fieldmap time values).  default: 'cubic'.
% <mcmask> (optional) is {MN SD} with the MN and SD inputs to makegaussian3d.m.
%   if supplied, we use these parameters to construct a binary 3D ellipse mask
%   that is used to determine which voxels contribute to the motion parameter
%   estimation (see defineellipse3d.m).  default to [] which means do nothing special.
% <maskoutnans> (optional) is
%   0 means do nothing special
%   1 means to zero out all of the data for any voxel that has NaN at any 
%     point in the EPI runs.  this option is recommended for maximum safety.
%   2 means to zero out all of the data in a given EPI run for any voxel that has
%     NaN at any point in that EPI run
%   default: 1.
% <epiignoremcvol> (optional) is like in motioncorrectvolumes.m.  default is []
%   which means do nothing special.  note that indices in <epiignoremcvol> should
%   be in reference to the data after volumes are dropped according
%   to <numepiignore>.
% <dformat> (optional) is 'single' | 'double'.  if you supply 'single', we will
%   attempt to use single format at various points in processing in order to
%   reduce memory usage.  default: 'double'.
% <epismoothfwhm> is a 3-element vector with the desired FWHM of a Gaussian filter.
%   if supplied, we smooth the EPI volumes right after slice time correction.
%   Default is [] which means do nothing.
%   Special case is [P X Y Z] where P is the integer number of pads to use
%   and X Y Z is the desired simulated voxel size -- in this case we use
%   replication to pad each side of the volumes and then use the <mode>==1 
%   case of smoothvolumes.m to use ideal Fourier filtering to achieve smoothing.
% <wantpushalt> is either [] which means do nothing special or a path to a record.mat
%   file from a previous call.  in this case, mcmask must not be {}, and the 
%   fmriquality stuff is skipped.  the effect of <wantpushalt> is to skip some 
%   processing --- specifically, we load in the previous <sliceshifts>, <mparams>,
%   and <smoothfieldmaps> and use them as-is.
%
% here's the short version:
%   we process the EPI data by (1) dropping the first few volumes (e.g. to avoid
% initial magnetization effects), (2) performing slice time correction via sinc
% interpolation, (3) calculating the center-of-mass of each slice at each time
% point with respect to the phase-encode dimension, band-pass filtering the results,
% and then shifting each slice along the phase-encode dimension to correct for
% the apparent motion, (4) undistorting the volumes based on the fieldmaps (after 
% unwrapping and smoothing the fieldmaps), (5) estimating motion parameters
% based on the slice-shifted and undistorted volumes, and then (6) performing the 
% slice-shifting, undistortion, and motion correction in one interpolation step 
% from the original EPI volumes. we return the corrected EPI data in the variable 
% <epis> in int16 format.  we write out figures illustrating all aspects of the 
% processing to <figuredir>.
%
% here's the long version:
%  1. if you supply in-plane volumes, we write them out as figures for inspection.
%     the in-plane volumes are individually contrast-normalized.  we also write out
%     versions of the in-plane volumes that are matched to the field-of-view and
%     resolution of the EPI data.  these versions are also individually
%     contrast-normalized.
%  2. for each EPI run, we drop volumes according to <numepiignore>.
%  3. for each EPI run, we perform slice time correction according to <episliceorder>,
%     interpolating each slice to the time of the first slice.  to obtain new values,
%     we use sinc interpolation, replicating the first and last time points to handle
%     edge issues.  (in the case where <episliceorder> is a cell vector of length 2,
%     we use pchip interpolation and change the TR of the data.)
%  4. for each EPI run, we compute the temporal SNR.  this is performed by regressing
%     out a line from each voxel's time-series, computing the absolute value of the
%     difference between successive time points, computing the median of these absolute
%     differences, dividing the result by the mean of the original time-series, and then
%     multiplying by 100.  negative values (which result when the mean is negative) are
%     explicitly set to NaN.  we write out the temporal SNR as figures for inspection,
%     using MATLAB's jet colormap.  high values (red) are good and correspond to a
%     temporal SNR of 0%.  low values (blue) are bad, and correspond to a temporal SNR
%     of 5%.  (note that it would make some sense to take the reciprocal of the computed
%     metric such that the mean signal level is in the numerator, but we leave it as-is
%     since we believe having the median absolute difference in the numerator is simpler.)
%  5. we write out the first and last volumes of each EPI run to the directory "EPIoriginal".
%     (the aggregate set of first volumes are contrast-normalized as a whole; the aggregate 
%     set of last volumes are contrast-normalized as a whole.)  we also write out the first
%     30 volumes of the first EPI run to the directory "MOVIEoriginal".  (the set of volumes
%     are contrast-normalized as a whole.)
%  6. if fieldmaps are provided, we write out the fieldmaps, fieldmap brains, and a histogram
%     of fieldmap values as figures for inspection.  we also write out versions of fieldmap
%     brains that are matched to the field-of-view and resolution of the EPI data.  we also
%     write out successive differences of the fieldmaps (e.g. 2-1, 3-2, 4-3, etc.), accounting
%     for phase wraparound.  the range of the fieldmap figures is -N Hz to N Hz, where N is 
%     1/(<fieldmapdeltate>/1000)/2.  the range of the fieldmap difference figures is -50 Hz 
%     to 50 Hz.  the fieldmap brain figures are individually contrast-normalized.
%  7. we use FSL's prelude utility to unwrap each fieldmap.
%  8. we write out the unwrapped fieldmaps as figures for inspection (range same as before).
%     we also write out figures that indicate the weighted mean of the unwrapped fieldmaps
%     in each slice.  the weights are given by the fieldmap brains.
%  9. we use local linear regression to smooth the unwrapped fieldmaps.  we obtain values of
%     the smoothed, unwrapped fieldmaps only at the locations of the EPI voxels.
% 10. we write out the smoothed, unwrapped fieldmaps as figures for inspection (range same as before).
%     we also write out versions that are matched to the field-of-view and the resolution
%     of the original fieldmap data.  we also write out "ALT" versions of these that have
%     a range of -N/3 Hz to N/3 Hz, in order to enhance visibility of small differences.
% 11. we undistort the first and last volumes of each run and write these out as figures
%     for inspection to the directory "EPIundistort".
% 12. if <sliceshiftband> is supplied, we calculate the center-of-mass of each EPI slice
%     at every time point and for both in-plane dimensions.  (EPI values are positively-
%     rectified to ensure valid center-of-mass calculations.)  then, we band-pass filter 
%     each center-of-mass time-series according to <sliceshiftband>.  we will use the 
%     band-pass filtered values obtained for the phase-encode direction to unshift each slice.
% 13. we write out figures illustrating the center-of-mass calculations.
% 13b. as a preliminary step, if motion correction is desired and the <mcmask> input is supplied,
%      we write out a figure illustrating the mask to be used in motion parameter estimation.
% 14a. if motion correction is desired, then we shift each slice (if applicable), undistort 
%      each volume (if applicable), estimate motion parameters, and then resample the volumes 
%      accounting for slice-shifting, undistortion, and motion correction.  for this final
%      step, we use only a single interpolation of the original EPI volumes (this avoids errors 
%      that could accumulate with multiple interpolations).  the final resampling is done using 
%      cubic interpolation.
% 14b. if motion correction is not desired, then we just shift each slice (if applicable) and
%      undistort each volume (if applicable), doing these in a single resampling step.
% 14c. if slice-shifting, undistortion, and motion correction are all not desired, then we 
%      do nothing.
% 15. if some correction was performed, we write out the first and last volumes of each run
%     as figures for inspection to the directory "EPIfinal" and the first 30 volumes of the
%     first EPI run to the directory "MOVIEfinal".  for these inspections, we detect
%     which of the three spatial dimensions has the fewest voxels and we make that dimension
%     the slice dimension.  the reason for this special handling is to ensure that case 2 of
%     <targetres> has reasonable inspection figures.
% 16. finally, we calculate:
%      <meanvolrun> as a cell vector of the mean volume of each EPI run (converted to int16)
%      <meanvol> as the mean volume aggregating over all EPI runs (converted to int16)
%      <validvolrun> as cell vector of logical volumes indicating which voxels had no NaNs for each EPI run
%      <validvol> as a logical volume indicating which voxels had no NaNs in any EPI run
%      <additionalvol> as a cell vector with various quantities (see #19 below)
%     then, as a last step, we zero out voxel data according to <maskoutnans>.
% 17. just before we finish up, we write out figures that illustrate the spatial quality of the
%     final corrected version of the EPI data.  (this is skipped when <targetres> is the 
%     cell vector case.)  the figures are written to a directory called "fmriquality".
%     see fmriquality.m for details on what the figures mean.  (note that the meanvolMATCH.png
%     figure is matched to the resolution and field-of-view of the _first_ in-plane volume.)
% 18. we save all of the various variables created in the processing to record.mat, excluding 
%     some of the big inputs to this function (namely, <inplanes>, <fieldmaps>, <fieldmapbrains>)
%     and the variables <fieldmapunwraps>, <sliceshifts>, and <finalfieldmaps>.
% 19. we return as output:
%     <epis> as the final corrected version of the EPI data.  the format is
%       same as the <epis> input, except that the numerical precision is now int16.
%       (in the special phase angle case, the input and output format are different (see above).)
%     <finalepisize> as a 3-element vector indicating the voxel size of the <epis> output.
%     <validvol> as a logical volume indicating which voxels had no NaNs in any EPI run.
%     <meanvol> as the mean volume aggregating over all EPI runs.
%     <additionalvol> as a cell vector with {A B} where
%        A is a volume with the median absolute difference (computed for the first run). invalid voxels get 0.
%        B is a volume with the temporal SNR (computed for the first run). invalid voxels get NaN.
%        please see computetemporalsnr.m for details on these quantities.
%
% notes:
% - we have deliberately offloaded format-specific loading and saving of data
%   to code outside of preprocessfmri.m.  a sample script that calls preprocessfmri.m
%   is given in preprocessfmri_standardscript.m.  this script assumes a certain
%   data setup.  for new setups, make your own scripts!
% - the undistortion process might assign bright voxel values to regions where fieldmap measurements are 
%   noisy (e.g. outside the brain).  be careful.  in particular, the assignment of bogus values to these 
%   regions may confuse (somewhat) the motion parameter estimation.
% - the value assigned to an EPI voxel may be NaN at any given point in a time-series (due to movement to 
%   outside of the field-of-view).  furthermore, because the <epis> variable is returned as int16, 
%   it cannot support NaN values.  thus, NaNs will show up as 0s in the <epis> variable.
%   for maximum safety, we recommend that you set the <maskoutnans> variable to 1, so that
%   every voxel will either have a complete set of data or will have its data set to all zeros.
% - we encapsulate the commands in this function with "dbstop if error" so that we can
%   perform troubleshooting if crashes occur.
%
% notes on phase angle case:
% - the general idea for handling <epis> in the case of phase angles is to represent the
%   angles as equidistant complex numbers, operate on the real and imaginary parts 
%   separately, and then re-convert back to equidistant complex numbers.  we have
%   to make sure all computations operate in a way consistent with this.  this includes
%   the temporal interpolation step, the smoothing step, and the spatial interpolation
%   step.  our standard phase angle format is complex unit-length numbers that have been
%   multiplied by 10000 and then converted to int16 (to save on precious memory!).
%
% notes on <targetres>:
% - in the case that <targetres> is {[A B C] [D E F] 1 H}, this is tricky.
%   what we do is to crop each EPI run to the smallest 3D box that contains all 
%   non-NaN values in the very first corrected volume of the run.  then, after
%   we have processed all EPI runs, we calculate the maximum 3D box that is still
%   contained within each of the EPI runs' 3D boxes.  we then crop all the EPI runs 
%   to this maximum 3D box.  note that this is an aggressive strategy that is 
%   happy to throw away voxels if these voxels would have had NaN values in one 
%   of the EPI runs anyway.  in the record.mat file, we save the variable
%   'voloffset' which is a 3-element vector of non-negative integers indicating
%   the offset relative to the original volume. for example, [100 0 40] means that 
%   the first voxel is actually (101,1,41) of the original volume.
%
% code dependencies:
% - SPM (for SPM's motion correction routines)
% - FSL (for FSL's prelude phase-unwrapping utility)
% - NIFTI_20110215 (for loading and saving NIFTI files)
% - ba_interp3 (for cubic interpolation of the unwrapped, smoothed fieldmaps and EPI volumes)
%
% note that .zip files containing the NIFTI_20110215 and ba_interp3 toolboxes are
% included in the kendrick repository for your convenience.
%
% technical notes (system administration stuff):
% - you may need to do something like
%     setenv MATLAB_SHELL /bin/tcsh
%   in your shell resource file in order for your environment variables to 
%   take effect when MATLAB calls the shell.  to check that things are working,
%   you can try typing
%     unix('prelude')
%   at the MATLAB prompt and see whether it can call prelude successfully.
% 
% history:
% 2017/11/29 - implement the [P X Y Z] case of <epismoothfwhm>
% 2016/12/27 - switch back to pchip temporal interpolation!
% 2016/08/09 - switch to spline temporal interpolation instead of cubic!!
% 2016/05/02 - add support for <wantpushalt> and the phase-angle case of <epis>; also,
%              the <additionalvol> stuff is now computed only for first run
% 2016/04/17 - add <additionalvol> output
% 2016/02/25 - expand flexibility of <episliceorder>
% 2016/02/05 - add <epismoothfwhm> input
% 2016/02/05 - the cell2 case of <episliceorder> now uses REPLICATION for the first and
%              last data points.  this changes previous behavior!!
% 2015/11/15 - implement the cell2 case of <episliceorder> and a few bug fixes
% 2015/02/28 - fix bug relating to zero-filling (would have crashed)
% 2014/11/26 - allow <episliceorder> to be the {X} case
% 2014/04/30 - allow <extratrans> to be the {X} case
% 2013/06/02 - back out previous change. (was buggy.).  must use NIFTI_20110215 apparently!
% 2013/05/28 - tweak to ensure compatibility with newer versions of the NIFTI toolbox
% 2011/08/07 - fix bug related to no inplanes being supplied (would have crashed)
% 2011/07/30 - add more option for episliceorder
% 2011/07/26 - after local regression of fieldmaps, force NaNs to be 0.
% 2011/07/26 - make <numepiignore> more flexible
% 2011/07/26 - allow epitr to be a vector
% 2011/07/15 - add <dformat>
% 2011/06/25 - add dbstop if error
% 2011/06/24 - write out some new figures ("ALT" fieldmaps)
% 2011/04/20 - fix bugs (it would have crashed); we now do save AFTER fmriquality.
% 2011/04/18 - write out fieldmap difference images
% 2011/04/17 - use the inplane input of fmriquality.m
% 2011/04/13 - implement <epiignoremcvol>; implement NaN case for <fmriqualityparams>
% 2011/04/13 - special [] case for <episliceorder>
% 2011/04/04 - fix minor bug (would have crashed)
% 2011/04/03 - implement input <maskoutnans>. note that the default is 1, so this changes previous behavior!
% 2011/03/29 - offload computetemporalsnr.m.
% 2011/03/29 - add writing of MOVIEoriginal and MOVIEfinal
% 2011/03/29 - add percentiles to fieldmap histograms
% 2011/03/28 - add many calls to reportmemoryandtime to aid in troubleshooting
% 2011/03/26 - implement <mcmask>.
% 2011/03/24 - implement <epiinplanematrixsize>.  the input format to this function changed!
% 2011/03/20 - meanvolrun and meanvol now returned as int16 (in the record.mat file). 
%              note that the entries in these variables that should be NaN will show up as 0.
%              also, remove validvxsrun and validvxs (no longer saved to record.mat file).
% 2011/03/19 - here's a big summary of all the recent changes.  hopefully things will be stable after this!:
%   big general changes: now we have the fieldmaptimes; add the inputs <sliceshiftband>, <fmriqualityparams>,
%     <fieldmaptimeinterp>; implement new fieldmap time interpolation capability; implement slice shifting; 
%     implement fmriquality figures.
%   for preprocessfmri.m, changes include: implement flattening option for <targetres>; exclude additional 
%     variables from being saved to the record.mat file; <epis> is now int16 format upon output; now we are 
%     much smarter handling of int16 data tricky issues; we now use ba_interp3 for interpolating the 
%     fieldmaps; new output variables meanvol and validvol; new input variables sliceshiftband,
%     fmriqualityparams,fieldmaptimeinterp; enforce single format for the finalfieldmaps and the sliceshifts
%     (in order to save memory)
%   for preprocessfmri_standard.m, we automatically make output directory if necessary; we save two new 
%     files (valid and mean).
%   for preprocessfmri_standardscript.m, we now check for matlabpool first before trying to open the pool, 
%     the default for episliceorder is now sequential, and we save new files valid.nii and mean.nii.  
%   for undistortvolumes.m, the changes were: switch to ba_interp3; use ba_interp3 for forward distortion 
%     (instead of MATLAB's interpolation); new output validvol; let's be explicit on data format issues; 
%     more flexible input; int16 for the output.
% 2011/03/19 - implement new special time-interpolation of fieldmaps
% 2011/03/18 - implement fmriquality
% 2011/03/18 - no longer save fieldmapunwraps
% 2011/03/16 - implement <sliceshiftband>
% 2011/03/15 - implement targetres H parameter
% 2011/03/14 - return epis as int16.
% 2011/03/10 - implement new case for <targetres>
% 2010/03/06 - first version!

% TODO:
% - impelement MI against inplane inspection?
% - can we make the script-level more automated?
% - kill all invalid voxels up front once and for all?
% - do we need to increase robustness of unwrap?
% - am i too aggressive on the edge interpolation strategy? (making edge voxels zero)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL CONSTANTS

% internal constants
tsnrmx = 5;          % max temporal SNR percentage (used in determining the color range)
numinchunk = 30;     % max images in chunk for movie
fmapdiffrng = [-50 50];  % range for fieldmap difference volumes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DBSTOP

dbstop if error;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREP

% input defaults
if ~exist('numepiignore','var') || isempty(numepiignore)
  numepiignore = 0;
end
if ~exist('motionreference','var') || isempty(motionreference)
  motionreference = [1 1];
end
if ~exist('motioncutoff','var') || isempty(motioncutoff)
  motioncutoff = 1/90;
end
if ~exist('extratrans','var') || isempty(extratrans)
  extratrans = eye(4);
end
if ~exist('targetres','var') || isempty(targetres)
  targetres = [];
end
if ~exist('sliceshiftband','var') || isempty(sliceshiftband)
  sliceshiftband = [];
end
if ~exist('fmriqualityparams','var') || isempty(fmriqualityparams)
  fmriqualityparams = {[] [] []};
end
if ~exist('fieldmaptimeinterp','var') || isempty(fieldmaptimeinterp)
  fieldmaptimeinterp = 'cubic';
end
if ~exist('mcmask','var') || isempty(mcmask)
  mcmask = [];
end
if ~exist('maskoutnans','var') || isempty(maskoutnans)
  maskoutnans = 1;
end
if ~exist('epiignoremcvol','var') || isempty(epiignoremcvol)
  epiignoremcvol = [];
end
if ~exist('dformat','var') || isempty(dformat)
  dformat = 'double';
end
if ~exist('epismoothfwhm','var') || isempty(epismoothfwhm)
  epismoothfwhm = [];
end
if ~exist('wantpushalt','var') || isempty(wantpushalt)
  wantpushalt = [];
end

% make cell if necessary
if ~iscell(inplanes)
  inplanes = {inplanes};
end
if ~iscell(inplanesizes)
  inplanesizes = {inplanesizes};
end
if ~iscell(fieldmaps)
  fieldmaps = {fieldmaps};  % now, it is either {vol}, {vol1 vol2 vol3...}, {A B}, or {}
end
if ~iscell(fieldmapbrains)
  fieldmapbrains = {fieldmapbrains};
end
if ~iscell(fieldmapsizes)
  fieldmapsizes = {fieldmapsizes};
end
if ~iscell(fieldmapsmoothing)
  fieldmapsmoothing = {fieldmapsmoothing};
end
if ~isempty(fieldmaps)  % disregard the {} case...
  if ~iscell(fieldmaps{1})
    fieldmaps = {fieldmaps};  % now, everything is {A B}, but B could be omitted or []
  end
  if length(fieldmaps) < 2
    fieldmaps{2} = [];  % insert B if necessary
  end
  if isempty(fieldmaps{2})
    fieldmaps{2} = 1:length(fieldmaps{1});  % handle [] case for B
  end
  fieldmaptimes = fieldmaps{2};  % split off B part into <fieldmaptimes>
  fieldmaps = fieldmaps{1};      % <fieldmaps> is now just the A part
  % ok, now, it is the case that either:
  % (1) fieldmaps is {} (which is the no undistortion case) and fieldmaptimes is not defined
  % (2) fieldmaps is A (cell vector of fieldmaps) and fieldmaptimes is fully specified
end
if ~iscell(epis)
  episissingle = 1;
  epis = {epis};
end

% calc
wantundistort = ~isempty(fieldmaps);
wantsliceshift = ~isempty(sliceshiftband);

% automatically repeat
if length(inplanesizes)==1
  inplanesizes = repmat(inplanesizes,[1 length(inplanes)]);
end
if length(fieldmapsizes)==1
  fieldmapsizes = repmat(fieldmapsizes,[1 length(fieldmaps)]);
end
if length(fieldmapdeltate)==1
  fieldmapdeltate = repmat(fieldmapdeltate,[1 length(fieldmaps)]);
end
if length(fieldmapunwrap)==1
  fieldmapunwrap = repmat({fieldmapunwrap},[1 length(fieldmaps)]);
end
if length(fieldmapsmoothing)==1
  fieldmapsmoothing = repmat(fieldmapsmoothing,[1 length(fieldmaps)]);
end
if length(epitr)==1
  epitr = repmat(epitr,[1 length(epis)]);
end
if (wantundistort || wantsliceshift) && length(epiphasedir)==1
  epiphasedir = repmat(epiphasedir,[1 length(epis)]);
end
if length(numepiignore)==1
  numepiignore = repmat(numepiignore,[1 length(epis)]);
end

% convert to special format
if ~iscell(numepiignore)
  numepiignore = cellfun(@(x) [x 0],num2cell(numepiignore),'UniformOutput',0);
end

% deal with more defaults
if isempty(epiinplanematrixsize)
  epiinplanematrixsize = sizefull(epis{1},2);
end

% deal with special extratrans case
if iscell(extratrans)
  dimdata = 1;
  dimtime = 2;
else
  dimdata = 3;
  dimtime = 4;
end

% calc
wantfigs = ~isempty(figuredir);
wantmotioncorrect = ~isequalwithequalnans(motionreference,NaN);
epidim = sizefull(epis{1},3);     % e.g. [64 64 20]
epifov = epidim .* episize;       % e.g. [128 128 40]

% convert fieldmapunwrap
for p=1:length(fieldmapunwrap)
  switch fieldmapunwrap{p}
  case 1
    fieldmapunwrap{p} = '-s -t 0';
  end
end

% convert sliceorder word cases
if ischar(episliceorder)
  switch episliceorder
  case 'sequential'
    episliceorder = 1:epidim;
  case 'interleaved'
    episliceorder = [1:2:epidim 2:2:epidim];
  case 'interleavedalt'
    if mod(epidim,2)==0
      episliceorder = [2:2:epidim 1:2:epidim];
    else
      episliceorder = [1:2:epidim 2:2:epidim];
    end
  otherwise
    error;
  end
end

% deal with defaults for episliceorder
if ~isempty(episliceorder) && iscell(episliceorder) && length(episliceorder)>=2
  if length(episliceorder{2})==1
    episliceorder{2} = repmat(episliceorder{2},[1 length(epis)]);
  end
  if length(episliceorder)<3
    episliceorder{3} = 0;
  end
  if length(episliceorder{3})==1
    episliceorder{3} = repmat(episliceorder{3},[1 length(epis)]);
  end
end

% deal with massaging epifieldmapasst [after this, epifieldmapasst will either be NaN or a fully specified cell vector]
if wantundistort
  if isempty(epifieldmapasst)
    if length(fieldmaps)==1  % if one fieldmap, use it for all
      epifieldmapasst = ones(1,length(epis));
    elseif length(fieldmaps)==length(epis)  % if equal number, assign 1-to-1
      epifieldmapasst = 1:length(fieldmaps);
    elseif length(fieldmaps)==length(epis)+1  % if one more fieldmap, then interpolate between successive
      epifieldmapasst = splitmatrix(flatten([1:length(fieldmaps)-1; 2:length(fieldmaps)]),2,2*ones(1,length(epis)));
    else
      error('<epifieldmapasst> cannot be [] when the number of fieldmaps is not one NOR the same as the number of EPI runs NOR the same as the number of EPI runs plus one');
    end
  end
  if ~iscell(epifieldmapasst) && ~isequalwithequalnans(epifieldmapasst,NaN)
    epifieldmapasst = num2cell(epifieldmapasst);
  end
  assert(isequalwithequalnans(epifieldmapasst,NaN) || (length(epifieldmapasst)==length(epis)));
end

% deal with targetres
if iscell(extratrans)
  targetres0 = [];    % targetres0 is used in specific function calls. this is a bit ugly.
else
  if isempty(targetres)
    targetres = epidim;
  end
  targetres0 = targetres(1:3);
end

% calc some special phase angle stuff
if isreal(epis{1})
  prefun = @(x) double(x);                            % regular case gets 1% normalization and gray colormap
  prerng = [];
  precmap = [];
else
  prefun = @(x) double(mod(angle(single(x)),2*pi));   % phase angle case gets fixed normalization and hsv colormap
  prerng = [0 2*pi*(255/256)];
  precmap = hsv(256);
end

% make figure dir
if wantfigs
  mkdirquiet(figuredir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO IT

  reportmemoryandtime;

% write out inplane volumes
if wantfigs
  fprintf('writing out inplane volumes for inspection...');
  for p=1:length(inplanes)
    imwrite(uint8(255*makeimagestack(inplanes{p},1)),sprintf('%s/inplane%02d.png',figuredir,p));
    imwrite(uint8(255*makeimagestack( ...
      processmulti(@imresizedifferentfov,inplanes{p},inplanesizes{p}(1:2),epidim(1:2),episize(1:2)), ...
      1)),sprintf('%s/inplaneMATCH%02d.png',figuredir,p));
  end
  fprintf('done.\n');
end

  reportmemoryandtime;

% drop the first few EPI volumes
fprintf('dropping EPI volumes (if requested).\n');
epis = cellfun(@(x,y) x(:,:,:,y(1)+1:end-y(2)),epis,numepiignore,'UniformOutput',0);

  reportmemoryandtime;

% slice time correct [NOTE: we may have to do in a for loop to minimize memory usage]
if ~isempty(episliceorder)
  fprintf('correcting for differences in slice acquisition times...');
  if iscell(episliceorder)
    if length(episliceorder)==1
      epis = cellfun(@(x,y) sincshift(x,repmat(reshape((1-y)/max(y),1,1,[]),[size(x,1) size(x,2)]),4), ...
                     epis,repmat({episliceorder{1}},[1 length(epis)]),'UniformOutput',0);
    else
      % this is the special case where we are changing the TR [pchip interpolation with padding]
      for p=1:length(epis)
        epistemp = cast([],class(epis{p}));
        for q=1:size(epis{p},3)  % process each slice separately
          temp0 = tseriesinterp(single(epis{p}(:,:,q,:)),epitr(p),episliceorder{2}(p),4,[], ...
                                -(((1-episliceorder{1}(q))/max(episliceorder{1})) * epitr(p)) - episliceorder{3}(p), ...
                                1,'pchip');
          if ~isreal(temp0)  % in the phase angle case, we have to revert back to true angles
            temp0 = int16(ang2complex(angle(temp0))*10000);
          end
          epistemp(:,:,q,:) = temp0;
        end
        epis{p} = epistemp;
      end
      clear epistemp temp0;
      epitr = episliceorder{2};  % we have a new TR, so change this!
    end
  else
    epis = cellfun(@(x,y) sincshift(x,repmat(reshape((1-y)/max(y),1,1,[]),[size(x,1) size(x,2)]),4), ...
                   epis,repmat({calcposition(episliceorder,1:max(episliceorder))},[1 length(epis)]),'UniformOutput',0);
  end
  fprintf('done.\n');
end

  reportmemoryandtime;

% smooth volumes if desired
if ~isempty(epismoothfwhm)
  fprintf('smoothing volumes...');
  if length(epismoothfwhm)==3
    epis = smoothvolumes(epis,episize,epismoothfwhm);
  else
    for zz=1:length(epis)
      temp = padarray(epis{zz},repmat(epismoothfwhm(1),[1 3]),'replicate','both');
      temp = smoothvolumes(temp,episize,epismoothfwhm(2:4),1);
      epis{zz} = temp(epismoothfwhm(1)+1:end-epismoothfwhm(1), ...
                      epismoothfwhm(1)+1:end-epismoothfwhm(1), ...
                      epismoothfwhm(1)+1:end-epismoothfwhm(1),:);
    end
    clear temp;
  end
  fprintf('done.\n');
end
if ~isreal(epis{1})  % in the phase angle case, we have to revert back to true angles
  for p=1:length(epis)
    epis{p} = int16(ang2complex(angle(single(epis{p})))*10000);
  end
end

  reportmemoryandtime;

% compute temporal SNR
  % this is a cell vector of 3D volumes.  values are percentages representing the median frame-to-frame difference
  % in units of percent signal.  (if the mean intensity is negative, the percent signal doesn't make sense, so
  % we set the final result to NaN.)  [if not enough volumes, some warnings will be reported.]
fprintf('computing temporal SNR...');
if isreal(epis{1})
  temporalsnr = cellfun(@computetemporalsnr,epis,'UniformOutput',0);
else
  temporalsnr = [];
end
fprintf('done.\n');

  reportmemoryandtime;

% write out EPI inspections
if wantfigs
  fprintf('writing out various EPI inspections...');

  % first and last of each run
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,1)),epis,'UniformOutput',0)),  sprintf('%s/EPIoriginal/image%%04da',figuredir),[],prerng,[],precmap);
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,end)),epis,'UniformOutput',0)),sprintf('%s/EPIoriginal/image%%04db',figuredir),[],prerng,[],precmap);

  % movie of first run
  viewmovie(prefun(epis{1}(:,:,:,1:min(30,end))),sprintf('%s/MOVIEoriginal/image%%04d',figuredir),[],prerng,[],precmap);

  % temporal SNR for each run
  if ~isempty(temporalsnr)
    for p=1:length(temporalsnr)
      imwrite(uint8(255*makeimagestack(tsnrmx-temporalsnr{p},[0 tsnrmx])),jet(256),sprintf('%s/temporalsnr%02d.png',figuredir,p));
    end
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% calc fieldmap stuff
fmapsc = 1./(fieldmapdeltate/1000)/2;  % vector of values like 250 (meaning +/- 250 Hz)

% write out fieldmap inspections
if wantfigs && wantundistort && isempty(wantpushalt)
  fprintf('writing out various fieldmap inspections...');

  % write out fieldmaps, fieldmaps brains, and histogram of fieldmap
  for p=1:length(fieldmaps)
  
    % write out fieldmap
    imwrite(uint8(255*makeimagestack(fieldmaps{p}/pi*fmapsc(p),[-1 1]*fmapsc(p))),jet(256),sprintf('%s/fieldmap%02d.png',figuredir,p));

    % write out fieldmap diff
    if p ~= length(fieldmaps)
      imwrite(uint8(255*makeimagestack(circulardiff(fieldmaps{p+1},fieldmaps{p},2*pi)/pi*fmapsc(p), ...
        fmapdiffrng)),jet(256),sprintf('%s/fieldmapdiff%02d.png',figuredir,p));
    end
  
    % write out fieldmap brain
    imwrite(uint8(255*makeimagestack(fieldmapbrains{p},1)),gray(256),sprintf('%s/fieldmapbrain%02d.png',figuredir,p));
    
    % write out fieldmap brain cropped to EPI FOV
    imwrite(uint8(255*makeimagestack(processmulti(@imresizedifferentfov,fieldmapbrains{p},fieldmapsizes{p}(1:2), ...
      epidim(1:2),episize(1:2)),1)),gray(256),sprintf('%s/fieldmapbraincropped%02d.png',figuredir,p));

    % write out fieldmap histogram
    figureprep; hold on;
    vals = prctile(fieldmaps{p}(:)/pi*fmapsc(p),[25 75]);
    hist(fieldmaps{p}(:)/pi*fmapsc(p),100);
    straightline(vals,'v','r-');
    xlabel('Fieldmap value (Hz)'); ylabel('Frequency');
    title(sprintf('Histogram of fieldmap %d; 25th and 75th percentile are %.1f Hz and %.1f Hz',p,vals(1),vals(2)));
    figurewrite('fieldmaphistogram%02d',p,[],figuredir);
  
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% unwrap fieldmaps
fieldmapunwraps = {};
if wantundistort && isempty(wantpushalt)
  fprintf('unwrapping fieldmaps if requested...');
  parfor p=1:length(fieldmaps)
  
    if ~isequal(fieldmapunwrap{p},0)
  
      % get temporary filenames
      tmp1 = tempname; tmp2 = tempname;
      
            %       % make a complex fieldmap and save to tmp1
            %       save_untouch_nii(make_ana(fieldmapbrains{p} .* exp(j*fieldmaps{p}),fieldmapsizes{p},[],32),tmp1);
      % make a complex fieldmap and save to tmp1
      save_nii(make_nii(fieldmapbrains{p} .* exp(j*fieldmaps{p}),fieldmapsizes{p},[],32),tmp1);

      % use prelude to unwrap, saving to tmp2
      unix_wrapper(sprintf('prelude -c %s -o %s %s; gunzip %s.nii.gz',tmp1,tmp2,fieldmapunwrap{p},tmp2));
      
      % load in the unwrapped fieldmap
      temp = load_nii(sprintf('%s.nii',tmp2));  % OLD: temp = readFileNifti(tmp2);
      
      % HACK (WHY IS THIS NECESSARY?)
      temp.img = flipdim(temp.img,1);
      
      % convert from radians centered on 0 to actual Hz
      fieldmapunwraps{p} = double(temp.img)/pi*fmapsc(p);
    
    else
  
      % convert from [-pi,pi] to actual Hz
      fieldmapunwraps{p} = fieldmaps{p}/pi*fmapsc(p);
  
    end
  
  end
  fprintf('done.\n');
end

  reportmemoryandtime;

% write out inspections of the unwrapping and additional fieldmap inspections
if wantfigs && wantundistort && isempty(wantpushalt)
  fprintf('writing out inspections of the unwrapping and additional inspections...');
  
  % write inspections of unwraps
  for p=1:length(fieldmaps)
    imwrite(uint8(255*makeimagestack(fieldmapunwraps{p},[-1 1]*fmapsc(p))),jet(256),sprintf('%s/fieldmapunwrapped%02d.png',figuredir,p));
  end

  % write slice-mean inspections
    % this is fieldmaps x slice-mean with the (weighted) mean of each slice in the fieldmaps:
  fmapdcs = catcell(1,cellfun(@(x,y) sum(squish(x.*abs(y),2),1) ./ sum(squish(abs(y),2),1),fieldmapunwraps,fieldmapbrains,'UniformOutput',0));
  figureprep; hold all;
  set(gca,'ColorOrder',jet(length(fieldmaps)));
  h = plot(fmapdcs');
  legend(h,mat2cellstr(1:length(fieldmaps)),'Location','NorthEastOutside');
  xlabel('Slice number'); ylabel('Weighted mean fieldmap value (Hz)');
  title('Inspection of fieldmap slice means');
  figurewrite('fieldmapslicemean',[],[],figuredir);

  fprintf('done.\n');
end

  reportmemoryandtime;

% use local linear regression to smooth the fieldmaps
smoothfieldmaps = cell(1,length(fieldmapunwraps));
if wantundistort && ~isequalwithequalnans(epifieldmapasst,NaN) && isempty(wantpushalt)
  fprintf('smooth the fieldmaps...');
  for p=1:length(fieldmapunwraps)
    if isnan(fieldmapsmoothing{p})
      smoothfieldmaps{p} = processmulti(@imresizedifferentfov,fieldmapunwraps{p},fieldmapsizes{p}(1:2),epidim(1:2),episize(1:2));
    else
      fsz = sizefull(fieldmaps{p},3);
      [xx,yy,zz] = ndgrid(1:fsz(1),1:fsz(2),1:fsz(3));
      [xi,yi] = calcpositiondifferentfov(fsz(1:2),fieldmapsizes{p}(1:2),epidim(1:2),episize(1:2));
      [xxB,yyB,zzB] = ndgrid(yi,xi,1:fsz(3));
      smoothfieldmaps{p} = nanreplace(localregression3d(xx,yy,zz,fieldmapunwraps{p},xxB,yyB,zzB,[],[],fieldmapsmoothing{p} ./ fieldmapsizes{p},fieldmapbrains{p},1),0,3);
    end
  end
  fprintf('done.\n');
end

  reportmemoryandtime;

% write out smoothed fieldmap inspections
if wantfigs && wantundistort && isempty(wantpushalt)
  fprintf('writing out smoothed fieldmaps...');

  % write out fieldmap and fieldmap resampled to match the original fieldmap
  for p=1:length(smoothfieldmaps)
    if ~isempty(smoothfieldmaps{p})
      todo = {{1 ''} {1/3 'ALT'}};
      for qqq=1:length(todo)
        imwrite(uint8(255*makeimagestack(smoothfieldmaps{p},todo{qqq}{1}*[-1 1]*fmapsc(p))),jet(256),sprintf('%s/fieldmapsmoothed%s%02d.png',figuredir,todo{qqq}{2},p));
        imwrite(uint8(255*makeimagestack(processmulti(@imresizedifferentfov,smoothfieldmaps{p},episize(1:2), ...
          sizefull(fieldmaps{p},2),fieldmapsizes{p}(1:2)),todo{qqq}{1}*[-1 1]*fmapsc(p))),jet(256), ...
          sprintf('%s/fieldmapsmoothedbacksampled%s%02d.png',figuredir,todo{qqq}{2},p));
      end
    end
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% deal with epifieldmapasst
finalfieldmaps = cell(1,length(epis));  % we need this to exist in all epi cases
if wantundistort
  fprintf('deal with epi fieldmap assignment and time interpolation...');

  % if push-alternative-data case, we have to load in smoothfieldmaps
  if ~isempty(wantpushalt)
    load(wantpushalt,'smoothfieldmaps');     % JUST-IN-TIME LOADING
  end

  % calculate the final fieldmaps [we use single to save on memory]
  if ~isequalwithequalnans(epifieldmapasst,NaN)
    for p=1:length(epifieldmapasst)
      if epifieldmapasst{p} ~= 0
        fn = epifieldmapasst{p};
        
        % if scalar, just use as-is, resulting in X x Y x Z
        if isscalar(fn)
          finalfieldmaps{p} = single(smoothfieldmaps{fn});
        
        % if two-element vector, do the interpolation, resulting in X x Y x Z x T     [[OUCH. THIS DOUBLES THE MEMORY USAGE]]
        else
          finalfieldmaps{p} = single(permute(interp1(fieldmaptimes,permute(catcell(4,smoothfieldmaps),[4 1 2 3]), ...
                                                     linspace(fn(1),fn(2),size(epis{p},4)),fieldmaptimeinterp,'extrap'),[2 3 4 1]));
        end
        
      end
    end
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% write out EPI undistort inspections
if wantfigs && wantundistort
  fprintf('writing out inspections of what the undistortion is like...');

  % undistort the first and last volume [NOTE: temp is int16, or complex int16]
  temp = cellfun(@(x,y,z) undistortvolumes(x(:,:,:,[1 end]),episize, ...
                 y(:,:,:,[1 end])*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2)), ...
                 z,[]),epis,finalfieldmaps,num2cell(epiphasedir),'UniformOutput',0);

  % inspect first and last of each run
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,1)),temp,'UniformOutput',0)),sprintf('%s/EPIundistort/image%%04da',figuredir),[],prerng,[],precmap);
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,2)),temp,'UniformOutput',0)),sprintf('%s/EPIundistort/image%%04db',figuredir),[],prerng,[],precmap);

  fprintf('done.\n');
end

  reportmemoryandtime;

% calculate center-of-mass stuff
if wantsliceshift

  if isempty(wantpushalt)

    % calculate center-of-mass after rectifying the epis.  each element of the cell vector is 1 x 2 x slices x time.
    com = cellfun(@(x) centerofmass(posrect(x),[1 2],2),epis,'UniformOutput',0);

    % apply band-pass filter
    combandpass = cellfun(@(x) reshape(zeromean(tsfilter(squish(x,3), ...
      constructbutterfilter1D(size(x,4),size(x,4)*sliceshiftband),[1 0 -1]),2),1,2,[],size(x,4)),com,'UniformOutput',0);
  
    % prepare the pixelshifts argument to undistortvolumes.m.  each element of the cell vector is epidim(1) x epidim(2) x slices x time.
    % we use single to save on memory.
    sliceshifts = cellfun(@(x,y) single(repmat(x(1,y,:,:),epidim(1:2))),combandpass,num2cell(abs(epiphasedir)),'UniformOutput',0);
  
  else
  
    load(wantpushalt,'sliceshifts');     % JUST-IN-TIME LOADING
    
  end
    
end

  reportmemoryandtime;

% write out inspections of center-of-mass stuff
if wantfigs && wantsliceshift && isempty(wantpushalt)
  fprintf('writing out inspections of slice-shifting stuff...');

  % process each run
  for p=1:length(com)
    figureprep([100 100 1000 600]); 
    for q=1:2
    
      temp =  bsxfun(@plus,1:epidim(3),zeromean(squish(com{p}(1,q,:,:),3)',1));  % time x slice
      temp2 =                          squish(combandpass{p}(1,q,:,:),3)';  % time x slice
      temp3 = bsxfun(@plus,1:epidim(3),squish(combandpass{p}(1,q,:,:),3)');  % time x slice

      subplot(1,6,(q-1)*3+1); hold on; 
      axis([0 size(com{p},4)+1 0 epidim(3)+1]);
      straightline(1:epidim(3),'h','k-'); 
      plot(temp);
      xlabel('Volume number'); ylabel('Center-of-mass offset in matrix units');
      title('Center-of-mass (original)');
      set(gca,'YTick',1:epidim(3),'YDir','reverse');
   
      subplot(1,6,(q-1)*3+2); hold on;
      axis([0 size(com{p},4)+1 0 epidim(3)+1]);
      straightline(1:epidim(3),'h','k-'); 
      plot(temp-temp2);
      xlabel('Volume number'); ylabel('Center-of-mass offset in matrix units');
      title('Center-of-mass (remove bandpass)');
      set(gca,'YTick',1:epidim(3),'YDir','reverse');

      subplot(1,6,(q-1)*3+3); hold on;
      axis([0 size(com{p},4)+1 0 epidim(3)+1]);
      straightline(1:epidim(3),'h','k-'); 
      plot(temp3);
      xlabel('Volume number'); ylabel('Center-of-mass offset in matrix units');
      title('Center-of-mass (bandpass)');
      set(gca,'YTick',1:epidim(3),'YDir','reverse');

    end
    figurewrite('centerofmass%02d',p,[],figuredir);
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% calc mcmaskvol and write out inspections
if wantmotioncorrect
  fprintf('writing out inspections of mcmaskvol (if applicable)...');

  % calculate mcmaskvol
  if isempty(mcmask)
    mcmaskvol = [];
  else
    mcmaskvol = double(makegaussian3d(epidim,mcmask{:}) > 0.5);
  end

  % inspect it
  if wantfigs && ~isempty(mcmask)
    imwrite(uint8(255*makeimagestack(mcmaskvol,[0 1])),gray(256),sprintf('%s/mcmaskvol.png',figuredir));
  end

  fprintf('done.\n');
end

  reportmemoryandtime;

% if we are doing motion correction, then...
fprintf('performing motion correction (if requested) and undistortion (if requested)...');
if wantmotioncorrect

  % if we are in the usual case (not pushing alternative data), then we have to calculate mparams and refvol
  if isempty(wantpushalt)
  
    % NOTE: the following two things could be put together into a single step...

    % slice-shift temporarily [NOTE: epistemp is int16 but gets converted to double/single]
    if wantsliceshift
      [epistemp,d,validvoltemp] = cellfun(@(x,y,z) undistortvolumes(x, ...
                         episize,y,z,[]),epis,sliceshifts,num2cell(abs(epiphasedir)),'UniformOutput',0);
      % yuck..  we have to explicitly convert to double/single and then set nan voxels to NaN
      for p=1:length(epistemp)
        epistemp{p} = squish(cast(epistemp{p},dformat),3);
        epistemp{p}(find(~validvoltemp{p}),:) = NaN;
        epistemp{p} = reshape(epistemp{p},sizefull(epis{p},4));
      end
    else
      epistemp = epis;
    end

    % undistort temporarily [NOTE: epistemp is int16 but gets converted to double/single]
    if wantundistort
      [epistemp,d,validvoltemp] = cellfun(@(x,y,z) undistortvolumes(x,episize, ...
        y*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2)),z,[]),epistemp,finalfieldmaps,num2cell(epiphasedir),'UniformOutput',0);
      % yuck..  we have to explicitly convert to double/single and then set nan voxels to NaN
      for p=1:length(epistemp)
        epistemp{p} = squish(cast(epistemp{p},dformat),3);
        epistemp{p}(find(~validvoltemp{p}),:) = NaN;
        epistemp{p} = reshape(epistemp{p},sizefull(epis{p},4));
      end
    end

    % estimate motion parameters from the slice-shifted and undistorted
    [epistemp,mparams,refvol] = motioncorrectvolumes(epistemp,cellfun(@(x,y) [x y],repmat({episize},[1 length(epis)]),num2cell(epitr),'UniformOutput',0), ...
      figuredir,motionreference,motioncutoff,[],1,[],[],mcmaskvol,epiignoremcvol,dformat);
    clear epistemp;
            %[epistemp,homogenizemask] = homogenizevolumes(epistemp,[99 1/4 2 2]);  % [],1
  
  % if we are in the pushing alternative data case, let's just load in mparams
  else
    
    load(wantpushalt,'mparams');     % JUST-IN-TIME LOADING

  end
  
  % finally, resample once (dealing with extratrans and targetres) [NOTE: epis is int16, or complex int16]
  if wantundistort
    if wantsliceshift
      [epis,voloffset,validvolrun] = cellfun(@(x,y0,y,z,w) undistortvolumes(x, ...
                     episize,bsxfun(@plus,sign(z)*y0, ...
                                    y*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2))),z,w(2:end,:),extratrans,targetres0), ...
                     epis,sliceshifts,finalfieldmaps,num2cell(epiphasedir),mparams,'UniformOutput',0);
    else
      [epis,voloffset,validvolrun] = cellfun(@(x,y,z,w) undistortvolumes(x,episize, ...
                     y*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2)),z,w(2:end,:),extratrans,targetres0), ...
                     epis,finalfieldmaps,num2cell(epiphasedir),mparams,'UniformOutput',0);
    end
  else
    if wantsliceshift
      [epis,voloffset,validvolrun] = cellfun(@(x,y0,z,w) undistortvolumes(x, ...
                     episize,y0,z,w(2:end,:),extratrans,targetres0), ...
                     epis,sliceshifts,num2cell(abs(epiphasedir)),mparams,'UniformOutput',0);
    else
      [epis,voloffset,validvolrun] = cellfun(@(x,w) undistortvolumes(x, ...
                     episize,[],[],w(2:end,:),extratrans,targetres0), ...
                     epis,mparams,'UniformOutput',0);
    end
  end

% if we're not doing motion correction then...
else

  % just slice-shift, undistort, and resample (dealing with extratrans and targetres) [NOTE: epis is int16, or complex int16]
  if wantundistort
    if wantsliceshift
      [epis,voloffset,validvolrun] = cellfun(@(x,y0,y,z) undistortvolumes(x, ...
                     episize,bsxfun(@plus,sign(z)*y0, ...
                                    y*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2))),z,[],extratrans,targetres0), ...
                     epis,sliceshifts,finalfieldmaps,num2cell(epiphasedir),'UniformOutput',0);
    else
      [epis,voloffset,validvolrun] = cellfun(@(x,y,z) undistortvolumes(x,episize, ...
                     y*(epireadouttime/1000)*(epidim(abs(z))/epiinplanematrixsize(2)),z,[],extratrans,targetres0), ...
                     epis,finalfieldmaps,num2cell(epiphasedir),'UniformOutput',0);
    end
  else
    if wantsliceshift
      [epis,voloffset,validvolrun] = cellfun(@(x,y0,z) undistortvolumes(x, ...
                     episize,y0,z,[],extratrans,targetres0), ...
                     epis,sliceshifts,num2cell(abs(epiphasedir)),'UniformOutput',0);
    else
      [epis,voloffset,validvolrun] = cellfun(@(x) undistortvolumes(x, ...
                     episize,[],[],[],extratrans,targetres0), ...
                     epis,'UniformOutput',0);
    end
  end

end
fprintf('done.\n');

%%% NOTE: at this point, for the phase angle case, we would normally want to 
%%% revert back to true angles, e.g. int16(ang2complex(angle(temp0))*10000)
%%% however, all the usage below of <epis> does not need this, so it would be 
%%% just unnecessary computational time.  so let's skip it.

  reportmemoryandtime;

% in the special case, we have to make sure all the epis
% are cropped the same way.  this is made necessary by the
% fact that we call undistortvolumes separately for each epi run.
% what a pain.
if iscell(targetres) && targetres{3}==1

  % 3-element vector with the indices of the very first voxel (playing it aggressively)
  beginix = max(catcell(1,voloffset)+1,[],1);

  % 3-element vector with the indices of the very last voxel (playing it aggressively)
  endix = min(catcell(1,voloffset) + catcell(1,cellfun(@(x) sizefull(x,3),epis,'UniformOutput',0)),[],1);
  
  % crop the epis and validvolrun
  epis = cellfun(@(x,y) x(beginix(1) - (y(1)+1) + 1 : end - (y(1)+size(x,1) - endix(1)), ...
                          beginix(2) - (y(2)+1) + 1 : end - (y(2)+size(x,2) - endix(2)), ...
                          beginix(3) - (y(3)+1) + 1 : end - (y(3)+size(x,3) - endix(3)),:),epis,voloffset,'UniformOutput',0);
  assert(all(cellfun(@(x) isequal(sizefull(x,3),sizefull(epis{1},3)),epis)));  % check that all runs have same matrix dimensions
  validvolrun = cellfun(@(x,y) x(beginix(1) - (y(1)+1) + 1 : end - (y(1)+size(x,1) - endix(1)), ...
                                 beginix(2) - (y(2)+1) + 1 : end - (y(2)+size(x,2) - endix(2)), ...
                                 beginix(3) - (y(3)+1) + 1 : end - (y(3)+size(x,3) - endix(3))),validvolrun,voloffset,'UniformOutput',0);
  
  % figure out the final voloffset
  voloffset = beginix - 1;
  
  % report to stdout
  fprintf('note: voloffset is %s and the 3D volume dimensions are %s.\n',mat2str(voloffset),mat2str(sizefull(epis{1},3)));

% otherwise, we just have to deal with preparing the <voloffset> variable for final output
else
  
  assert(all(cellfun(@(x) isequal(x,[0 0 0]),voloffset)));
  voloffset = [0 0 0];
  
end

  reportmemoryandtime;

% write out EPI final inspections
if wantfigs && (wantmotioncorrect || wantundistort || wantsliceshift) && ~iscell(extratrans)
  fprintf('writing out inspections of final EPI results...');

  % inspect first and last of each run
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,1)),epis,'UniformOutput',0)),  sprintf('%s/EPIfinal/image%%04da',figuredir),[],prerng,1,precmap);
  viewmovie(catcell(4,cellfun(@(x) prefun(x(:,:,:,end)),epis,'UniformOutput',0)),sprintf('%s/EPIfinal/image%%04db',figuredir),[],prerng,1,precmap);

  % inspect movie of first run
  viewmovie(prefun(epis{1}(:,:,:,1:min(30,end))),sprintf('%s/MOVIEfinal/image%%04d',figuredir),[],prerng,1,precmap);

  fprintf('done.\n');
end

  reportmemoryandtime;

% final EPI calculations  [note: mean of int16 produces double! this is good, except for the NaN issue]
fprintf('performing final EPI calculations...');
if isreal(epis{1})
  meanvolrun = cellfun(@(x) int16(mean(x,dimtime)),epis,'UniformOutput',0);          % mean of each run
  meanvol = int16(mean(catcell(dimtime,epis),dimtime));                              % mean over all runs
else
  meanvolrun = [];
  meanvol = [];
end
validvolrun = validvolrun;                                                         % logical of which voxels have no nans (in each run)
validvol = all(catcell(dimtime,validvolrun),dimtime);                              % logical of which voxels have no nans (over all runs)
  % deal with temporal SNR
if isreal(epis{1})
  [tsnr,~,mad] = computetemporalsnr(single(epis{1}),dimtime);
  additionalvol = {mad tsnr};
else
  additionalvol = {[]  []  };
end
if iscell(extratrans)
  finalepisize = [];
else
  if iscell(targetres)
    finalepisize = targetres{2};
  else
    finalepisize = epifov ./ targetres;                                        % size in mm of a voxel in the final EPI version
  end
end
  % some final adjustments for invalid voxels
if ~isempty(meanvolrun)
  meanvolrun = cellfun(@(x,y) copymatrix(x,~y,0),meanvolrun,validvolrun,'UniformOutput',0);
end
if ~isempty(meanvol)
  meanvol(~validvol) = 0;                  % invalid gets mean 0
end
if ~isempty(additionalvol{1})
  additionalvol{1}(~validvol) = 0;         % invalid gets mad 0
end
if ~isempty(additionalvol{2})
  additionalvol{2}(~validvol) = NaN;       % invalid gets tsnr NaN
end
fprintf('done.\n');

  reportmemoryandtime;

% in the case of phase angles, we have to do some final post-processing of epis
if ~isreal(epis{1})
  for p=1:length(epis)
    epis{p} = mod(int16(mod(angle(single(epis{p})),2*pi) * (4095/(2*pi))),4095);  % reverse is double(x) * (2*pi/4095)
  end
end

% zero out data
fprintf('zeroing out data for bad voxels...');
switch maskoutnans
case 0
case 1
  epis = cellfun(@(x) copymatrix(x,repmat(~validvol,[ones(1,dimdata) size(x,dimtime)]),0),epis,'UniformOutput',0);
case 2
  epis = cellfun(@(x,y) copymatrix(x,repmat(~y,[ones(1,dimdata) size(x,dimtime)]),0),epis,validvolrun,'UniformOutput',0);
end
fprintf('done.\n');

  reportmemoryandtime;

% clear out some variables!!!
clear xx yy zz xxB yyB zzB temp;
clear fieldmaps fieldmapbrains fieldmapunwraps;
clear sliceshifts finalfieldmaps;
clear prefun
  
% do fMRI quality
if wantfigs && ~iscell(targetres) && ~iscell(extratrans) && ~isequalwithequalnans(fmriqualityparams,NaN) && isempty(wantpushalt)
  fprintf('calling fmriquality.m on the epis...');
  if ~isempty(inplanes)
    inplaneextra = {sizefull(inplanes{1},2) inplanesizes{1}(1:2)};
  else
    inplaneextra = [];
  end
  fmriquality(epis,episize,fullfile(figuredir,'fmriquality'),fmriqualityparams{:},inplaneextra);  % note that NaNs are present and may be in weird places...
  fprintf('done with fmriquality.m.\n');
end

  reportmemoryandtime;

% save record
if ~isempty(figuredir)
  fprintf('saving record.mat...');
  clear inplanes;
  saveexcept(fullfile(figuredir,'record.mat'),'epis');  % ignore this big variable, but we need it upon function completion
  fprintf('done.\n');
end

  reportmemoryandtime;

% prepare epis in the special flattening case
if iscell(targetres) && targetres{4}==1
  fprintf('preparing epis for special flattening...');
  for p=1:length(epis)
    epis{p} = reshape(epis{p},[prod(sizefull(epis{p},3)) 1 1 size(epis{p},4)]);
    epis{p} = epis{p}(find(validvol),:,:,:);
  end
  fprintf('done.\n');
end

  reportmemoryandtime;

% prepare output
if exist('episissingle','var') && episissingle
  epis = epis{1};
end

  reportmemoryandtime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DBSTOP

dbclear if error;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK:

% internal constants
%fmapseedfctr = 5;    % number of seeds along each dimension of the fieldmap
%fmapptile = 99;      % percentile of absolute value of fieldmap values for scaling purposes
%fmapbptile = 99;     % percentile of fieldmap brain values for scaling purposes
% fmapmaskpct = 99;    % percentile of fieldmap brain for the mask
% fmapmaskfctr = 0.1;  % multiplier on the percentile for the mask



  %     % try many starting seeds for the unwrapping
  %     ii = resamplingindices(1,size(fieldmaps{p},1),fmapseedfctr);
  %     jj = resamplingindices(1,size(fieldmaps{p},2),fmapseedfctr);
  %     kk = resamplingindices(1,size(fieldmaps{p},3),fmapseedfctr);
  %     fieldmapedgepower = zeros(fmapseedfctr,fmapseedfctr,fmapseedfctr);
  %     for i1=1:length(ii)
  %       for i2=1:length(jj)
  %         for i3=1:length(kk)
  %           temp = robustunwrap([ii(i1) jj(i2) kk(i3)],fieldmaps{p}+pi,fieldmapbrains{p});
  %           fieldmapedgepower(i1,i2,i3) = sum(flatten(detectedges(temp,1).^2 .* abs(fieldmapbrains{p})));  % abs for safety
  %         end
  %       end
  %     end
  %     
  %     % which was best
  %     [d,ix] = min(fieldmapedgepower(:));
  %     [i1,i2,i3] = ind2sub([fmapseedfctr fmapseedfctr fmapseedfctr],ix);
  %     
  %     % do the unwrapping one last time
  %     fieldmapunwraps{p} = robustunwrap([ii(i1) jj(i2) kk(i3)],fieldmaps{p}+pi,fieldmapbrains{p});
  %     
  %     % convert from [0,2*pi] to actual Hz
  %     fieldmapunwraps{p} = (fieldmapunwraps{p}-pi)/pi*fmapsc(p);

      % % calculate fieldmap-related ranges
      % if wantundistort
      %   fmapmx = prctile(abs(catcell(1,cellfun(@(x) x(:),fieldmaps,'UniformOutput',0))),fmapptile);    % [-A A] range for fieldmaps
      %   fmapbmx = prctile(catcell(1,cellfun(@(x) x(:),fieldmapbrains,'UniformOutput',0)),fmapbptile);  % [0 B] range for fieldmap brains
      % end
%[other strategies?: voting (median). only the best ones. std dev. zoom before selecting seeds.]

%   % inspect special volumes of each run
%   viewmovie(catcell(4,cellfun(@(x) double(x(1:max(1,round(size(x,1)/64)):end,1:max(1,round(size(x,2)/64)):end,1:max(1,round(size(x,3)/64)):end, ...
%     round(linspace(1,size(x,4),9)))),epis,'UniformOutput',0)), ...
%     sprintf('%s/EPIfinalalt/image%%04d',figuredir),[],[],1);
% 

% validvxsrun = cellfun(@(x) flatten(find(x)),validvolrun,'UniformOutput',0);  % vector of indices of valid voxels (in each run)
% validvxs = flatten(find(validvol));                                          % vector of indices of valid voxels (over all runs)
