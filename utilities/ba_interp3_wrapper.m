function f = ba_interp3_wrapper(vol,coords,interptype)

% function f = ba_interp3_wrapper(vol,coords,interptype)
%
% <vol> is a 3D matrix
% <coords> is 3 x N with the matrix coordinates to interpolate at.
%   one or more of the entries can be NaN.
% <interptype> (optional) is 'nearest' | 'linear' | 'cubic'.  default: 'cubic'.
%
% this is a convenient wrapper for ba_interp3.  the main problem with
% normal calls to ba_interp3 is that it assigns values to interpolation
% points that lie outside the original data range.  what we do is to 
% ensure that coordinates that are outside the original field-of-view
% (i.e. if the value along a dimension is less than 1 or greater than
% the number of voxels in the original volume along that dimension)
% are returned as NaN and coordinates that have any NaNs are returned
% as NaN.
%
% we automatically convert <vol> and <coords> to double before calling
% ba_interp3.m because it requires double format.  the output from
% this function is also double.
%
% history:
% 2011/03/19 - be explicit on double conversion and double output.
%
% example:
% vol = -getsamplebrain(1);
% [xx,yy,zz] = ndgrid(1:size(vol,1),1:size(vol,2),1:size(vol,3));
% xx = xx + 2;
% yy = yy + 8;
% newvol = reshape(ba_interp3_wrapper(vol,[flatten(xx); flatten(yy); flatten(zz)]),size(vol));
% figure; imagesc(makeimagestack(vol));
% figure; imagesc(makeimagestack(newvol));

% input
if ~exist('interptype','var') || isempty(interptype)
  interptype = 'cubic';
end

% bad locations must get set to NaN
bad = any(isnan(coords),1);
coords(:,bad) = 1;

% out of range must become NaN, too   [ba_interp3 replicates the edge values which is not what we want]
bad = bad | ...
  coords(1,:) < 1 | coords(1,:) > size(vol,1) | ...
  coords(2,:) < 1 | coords(2,:) > size(vol,2) | ...
  coords(3,:) < 1 | coords(3,:) > size(vol,3);

% resample the volume
f = copymatrix(ba_interp3(double(vol),double(coords(2,:)),double(coords(1,:)),double(coords(3,:)),interptype),bad,NaN);
