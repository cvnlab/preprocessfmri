% history:
% 2015/02/28 - circumvent special fields if fieldmapB0files is empty anyway
% 2013/04/14 - determine epiphasedir on a run-by-run basis (so as to allow epiphasedir to
%              be different on different runs)
% 2013/03/15 - add 60 seconds pause; fix bug; save mean.nii and valid.nii files now
% 2013/03/08 - move matlabpool to the script; special header variable names changed
% 2013/03/08 - add in text reporting of fieldmapB0files and fieldmapMAGfiles
% 2013/03/04 - add back numepiignore
% 2013/03/04 - automate epiinplanematrixsize, epiphasedir, epireadouttime based on CNI header information
% 2013/02/27 - first version

% this file is called by preprocessfmri_CNIscript.m

if isempty(fieldmapB0files)
  fprintf('no fieldmapB0files were specified (VERIFY THAT THIS IS CORRECT).\n\n');
else
  fprintf('the following are the fieldmapB0files that we found (VERIFY THAT THIS IS CORRECT):\n');
  cellfun(@(x) fprintf(['  ' x '\n']),fieldmapB0files);
  fprintf('\n');
end

if isempty(fieldmapMAGfiles)
  fprintf('no fieldmapMAGfiles were specified (VERIFY THAT THIS IS CORRECT).\n\n');
else
  fprintf('the following are the fieldmapMAGfiles that we found (VERIFY THAT THIS IS CORRECT):\n');
  cellfun(@(x) fprintf(['  ' x '\n']),fieldmapMAGfiles);
  fprintf('\n');
end

if isempty(inplanefilenames)
  fprintf('no inplanefilenames were specified (VERIFY THAT THIS IS CORRECT).\n\n');
else
  fprintf('the following are the inplanefilenames that we found (VERIFY THAT THIS IS CORRECT):\n');
  cellfun(@(x) fprintf(['  ' x '\n']),inplanefilenames);
  fprintf('\n');
end

if isempty(epifilenames)
  fprintf('no epifilenames were specified (VERIFY THAT THIS IS CORRECT).\n\n');
else
  fprintf('the following are the epifilenames that we found (VERIFY THAT THIS IS CORRECT):\n');
  cellfun(@(x) fprintf(['  ' x '\n']),epifilenames);
  fprintf('\n');
end

fprintf('***** Please verify that the above files are correct.  We will proceed in 60 seconds. *****\n\n');
pause(60);

reportmemoryandtime;

% load Inplane files
fprintf('loading inplane data...');
inplanes = {}; inplanesizes = {};
for p=1:length(inplanefilenames)
  ni = load_untouch_nii(gunziptemp(inplanefilenames{p}));
  inplanes{p} = double(ni.img);
  inplanesizes{p} = ni.hdr.dime.pixdim(2:4);
  clear ni;
end
if exist('inplanehackfun','var')  % HRM. HACKY.
  inplanes = cellfun(inplanehackfun,inplanes,'UniformOutput',0);
end
fprintf('done (loading inplane data).\n');

reportmemoryandtime;

if ~exist('dformat','var')
  dformat = [];
end

% interactive prompt for mcmask
if iscell(mcmask) && isempty(mcmask)
  fprintf('loading first EPI run so that we can define an ellipse...');
  ni = load_untouch_nii(gunziptemp(epifilenames{1}));
  tempepi = double(ni.img);
  fprintf('done (loading first EPI run).\n');
  [d,tempmn,tempsd] = defineellipse3d(tempepi(:,:,:,1),[],0);
  mcmask = {tempmn tempsd};
  fprintf('mcmask = %s;\n',cell2str(mcmask));
  clear ni tempepi;
end

reportmemoryandtime;

% load EPI files
fprintf('loading EPI data...');
epis = {}; episizes = {}; epitr = {};
epiphasedir = [];
for p=1:length(epifilenames)
  ni = load_untouch_nii(gunziptemp(epifilenames{p}));
  epis{p} = single(ni.img);
  episizes{p} = ni.hdr.dime.pixdim(2:4);
  epitr{p} = ni.hdr.dime.pixdim(5);
  
  if ~isempty(fieldmapB0files)

    % if this is the first EPI file, then attempt to learn some information based on
    % the information in the header.
    if p==1

      % CNI puts some helpful information in the "descrip" field.  here,
      % we simply evaluate it (since it is valid MATLAB code).  the variables
      % defined are as follows:
      %   ec is the echo spacing (read-out time per PE line) in milliseconds
      %   rp is the ASSET/ARC acceleration factor in the phase-encode dimension
      %   acq is the acquisition matrix size (freq x phase)
      eval(ni.hdr.hist.descrip);  % this should define [ec, rp, acq]

      % this should be [A B] where A and B are the in-plane frequency-encode
      % and phase-encode matrix sizes, respectively.  can be [] in which case 
      % we default to the size of the first two dimensions of the EPI data.
      epiinplanematrixsize = acq;
      fprintf('*** epiinplanematrixsize determined to be %s.\n',mat2str(epiinplanematrixsize));
    
      % what is the total readout time in milliseconds for an EPI slice?
      % (note that 'Private_0043_102c' in the dicominfo of the EPI files gives the time per phase-encode line in microseconds.
      % I confirmed that this time is correct by checking against the waveforms displayed by plotter.)
      epireadouttime = ec*acq(2)/rp;  % divide by 2 if you are using 2x acceleration
      fprintf('*** epireadouttime determined to be %.5f ms * %d lines / %d acceleration.\n',ec,acq(2),rp);
    
    end

    % what is the phase-encode direction for the EPI runs? (see preprocessfmri.m for details.)
    % up-down in the images is 1 or -1 in our convention; left-right in the images is 2 or -2 
    % in our convention.  you should always check the sanity of the results!
    % NOTE: this attempts to learn this information from the NIFTI.
    %       if you ever flip the phase-encode direction, you will need to multiply
    %       the following by -1.
    epiphasedir(p) = bitand(uint8(3),bitshift(uint8(ni.hdr.hk.dim_info),-2));
    fprintf('*** epiphasedir for run %d determined to be %d.\n',p,epiphasedir(p));

  else
    epiinplanematrixsize = [];
    epireadouttime = [];
    epiphasedir = [];
  end

  clear ni;
end
fprintf('done (loading EPI data).\n');

reportmemoryandtime;

% load fieldmap data
fprintf('loading fieldmap data...');
fieldmaps = {}; fieldmapsizes = {}; fieldmapbrains = {};
for p=1:length(fieldmapB0files)
  ni = load_untouch_nii(gunziptemp(fieldmapB0files{p}));
  fieldmaps{p} = double(ni.img) * pi / (1/(fieldmapdeltate/1000)/2) ;  % convert to range [-pi,pi]
  fieldmapsizes{p} = ni.hdr.dime.pixdim(2:4);
  ni = load_untouch_nii(gunziptemp(fieldmapMAGfiles{p}));
  fieldmapbrains{p} = double(ni.img(:,:,:,1));  % JUST USE FIRST VOLUME
  clear ni;
end
fprintf('done (loading fieldmap data).\n');

reportmemoryandtime;

% deal with upsampling
fprintf('resampling fieldmap data if necessary...');
  % defaults for backwards-compatibility:
  if ~exist('fieldmapslicefactor','var')
    fieldmapslicefactor = [];
  end
if ~isempty(fieldmapslicefactor)
  if length(fieldmapsizes)==1
    fieldmapsizes = repmat(fieldmapsizes,[1 length(fieldmaps)]);  % make full just to make life easier
  end
  for p=1:length(fieldmaps)
    fieldmaps{p} = upsamplematrix(fieldmaps{p},[1 1 fieldmapslicefactor],[],[],'nearest');
    fieldmapbrains{p} = upsamplematrix(fieldmapbrains{p},[1 1 fieldmapslicefactor],[],[],'nearest');
    fieldmapsizes{p}(3) = fieldmapsizes{p}(3) / fieldmapslicefactor;
  end
end
fprintf('done (resampling fieldmap data).\n');

reportmemoryandtime;

% start parallel MATLAB to speed up execution.
if matlabpool('size')==0
  matlabpool open;
end

% do the pre-processing
  % defaults for backwards-compatibility:
  if ~exist('maskoutnans','var')
    maskoutnans = [];
  end
  if ~exist('epiignoremcvol','var')
    epiignoremcvol = [];
  end
fprintf('calling preprocessfmri...');
[epis,finalepisize,validvol,meanvol] = preprocessfmri(figuredir,inplanes,inplanesizes, ...
  {fieldmaps fieldmaptimes},fieldmapbrains,fieldmapsizes,fieldmapdeltate,fieldmapunwrap,fieldmapsmoothing, ...
  epis,episizes{1},epiinplanematrixsize,cell2mat(epitr),episliceorder, ...
  epiphasedir,epireadouttime,epifieldmapasst, ...
  numepiignore,motionreference,motioncutoff,extratrans,targetres, ...
  sliceshiftband,fmriqualityparams,fieldmaptimeinterp,mcmask,maskoutnans,epiignoremcvol,dformat);
fprintf('done (calling preprocessfmri).\n');

reportmemoryandtime;

% save it
fprintf('saving data...');
mkdirquiet(stripfile(savefile));
for p=1:length(epis)
  if iscell(targetres) && length(targetres) >= 4 && targetres{4}==1
    fprintf('for EPI run %d, we have %d time points and %d valid voxels.\n',p,size(epis{p},4),size(epis{p},1));
    savebinary(sprintf(savefile,p),'int16',squish(int16(epis{p}),3)');  % special flattened format: time x voxels
  else
    ni = load_untouch_nii(gunziptemp(epifilenames{p}));
    assert(isequal(sizefull(ni.img,3),sizefull(epis{p},3)));
    ni.img = cast(epis{p},class(ni.img));
    ni.hdr.dime.dim(5) = size(ni.img,4);  % since the number of volumes may have changed
    save_untouch_nii(ni,sprintf(savefile,p));
    
    % save special files
    if p==1
      ni.img = cast(validvol,class(ni.img));
      ni.hdr.dime.dim(5) = 1;
      save_untouch_nii(ni,sprintf([stripfile(savefile) '/valid.nii']));

      ni.img = cast(meanvol,class(ni.img));
      ni.hdr.dime.dim(5) = 1;
      save_untouch_nii(ni,sprintf([stripfile(savefile) '/mean.nii']));
    end

    clear ni;
  end
end
fprintf('done (saving data).\n');

reportmemoryandtime;
