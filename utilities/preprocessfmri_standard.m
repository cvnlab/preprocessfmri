% history:
% 2011/08/07 - fix bugs related to no inplanes and/or no fieldmaps (would have crashed)
% 2011/07/30 - load first EPI run quickly for mcmask definition
% 2011/07/28 - mcmask report is an actual matlab command now
% 2011/07/27 - allow fieldmapfiles to be NaN (user specification)
% 2011/07/27 - now load epitr directly from the dicoms!!
% 2011/07/26 - add fieldmapslicefactor
% 2011/07/26 - epitr more general
% 2011/07/26 - automatically take only the first N volumes of fieldmapbrains
% 2011/07/21 - add support for Siemens data
% 2011/07/21 - remove i%06d.dcm hint
% 2011/07/15 - add <dformat>
% 2011/06/24 - move EPI loading before fieldmap so that you can define the ellipse interactively quickly
% 2011/04/13 - add new input <epiignoremcvol>
% 2011/04/11 - allow fieldmapfiles to be the cached .mat files.
% 2011/04/04 - previously, the TR being saved in the .nii files was 1.  we now save
%              the actual TR value (as given by <epitr>) into the .nii files.

% this file is called by preprocessfmri_standardscript.m

if isempty(fieldmapfiles)
  fprintf('no fieldmapfiles were specified (VERIFY THAT THIS IS CORRECT).\n\n');
elseif iscell(fieldmapfiles)
  fprintf('the following are the fieldmapfiles that we found (VERIFY THAT THIS IS CORRECT):\n');
  if iscell(fieldmapfiles{1})
    cellfun(@(x) fprintf(['  fieldmap: ' x{1} ', magnitude: ' x{2} '\n']),fieldmapfiles);
  else
    cellfun(@(x) fprintf(['  ' x '\n']),fieldmapfiles);
  end
  fprintf('\n');
elseif isequalwithequalnans(fieldmapfiles,NaN)
  fprintf('fieldmapfiles was set to NaN, indicating that the user specified the data.\n');
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

reportmemoryandtime;

% load Inplane DICOM directories
fprintf('loading inplane data...');
[inplanes,inplanesizes,inplaneinplanematrixsizes] = dicomloaddir(inplanefilenames);  %,'i%06d.dcm');
if exist('inplanehackfun','var')  % HRM. HACKY.
  inplanes = cellfun(inplanehackfun,inplanes,'UniformOutput',0);
end
fprintf('done (loading inplane data).\n');

reportmemoryandtime;

if ~exist('epinumonly','var')
  epinumonly = [];
end
if ~exist('epidesiredinplanesize','var')
  epidesiredinplanesize = [];
end
if ~exist('epiphasemode','var')
  epiphasemode = [];
end
if ~exist('dformat','var')
  dformat = [];
end

% interactive prompt for mcmask
if iscell(mcmask) && isempty(mcmask)
  fprintf('loading first EPI run so that we can define an ellipse...');
  [tempepi] = dicomloaddir(epifilenames(1),[],epinumonly,epidesiredinplanesize,epiphasemode,dformat);
  fprintf('done (loading first EPI run).\n');
  [d,tempmn,tempsd] = defineellipse3d(tempepi{1}(:,:,:,1),[],0);
  mcmask = {tempmn tempsd};
  fprintf('mcmask = %s;\n',cell2str(mcmask));
  clear tempepi;
end

reportmemoryandtime;

% load EPI DICOM directories
fprintf('loading EPI data...');
[epis,episizes,epiinplanematrixsizes,epitr] = dicomloaddir(epifilenames,[],epinumonly,epidesiredinplanesize,epiphasemode,dformat);
%         epis = {}; episizes = {};
%         for p=1:length(epifilenames)
%           temp = readFileNifti(epifilenames{p});
%           epis{p} = double(temp.data);
%           episizes{p} = temp.pixdim(1:3);
%           clear temp;
%         end
%         epiindex = feval(epiindexfun,epis);
fprintf('done (loading EPI data).\n');

reportmemoryandtime;

% load fieldmap data
fprintf('loading fieldmap data...');

%%% this is the special user-specification case
if isequalwithequalnans(fieldmapfiles,NaN)
  assert(exist('fieldmaps','var') & exist('fieldmapsizes','var') & exist('fieldmapbrains','var'));  

%%% this is the DICOM case (first DICOM is the fieldmap; second DICOM is the magnitude brain)
elseif ~isempty(fieldmapfiles) && iscell(fieldmapfiles{1})
  [fieldmaps,fieldmapsizes] = dicomloaddir(cellfun(@(x) x{1},fieldmapfiles,'UniformOutput',0));  % load fieldmaps
  fieldmaps = cellfun(@(x) normalizerange(x,-pi,pi),fieldmaps,'UniformOutput',0);  % normalize min to -pi and max to pi
  [fieldmapbrains] = dicomloaddir(cellfun(@(x) x{2},fieldmapfiles,'UniformOutput',0));  % load magnitude brains
  fieldmapbrains = cellfun(@(x,y) x(:,:,1:size(y,3)),fieldmapbrains,fieldmaps,'UniformOutput',0);  % use first set

%%% this is the raw spiral *.7 case
else

  % reconstruct and prepare spiral fieldmaps [note that we save results to .mat files and attempt to re-use the results if they exist]
    % NOTE: making this parallel is non-trivial it seems due to weird errors.  but it seems fast enough.
  fieldmaps = {};
  fieldmapbrains = {};
  for p=1:length(fieldmapfiles)
  
    % if it is a .mat file, load from it
    if isequal(getextension(fieldmapfiles{p}),'.mat')
      load(fieldmapfiles{p},'fieldmap','fieldmapbrain');
  
    % if a .mat file exists, load from it
    elseif exist([fieldmapfiles{p} '.mat'],'file')
      load([fieldmapfiles{p} '.mat'],'fieldmap','fieldmapbrain');
  
    % otherwise, we have to recon it
    else
    
      % recon
      [fieldmap,fieldmapbrain] = reconstructspiralfieldmap(fieldmapfiles{p});
      fieldmap = fieldmap / (1/(fieldmapdeltate/1000)/2) * pi;  % convert to range [-pi,pi]
      
      % save results
      save([fieldmapfiles{p} '.mat'],'fieldmap','fieldmapbrain');
    
    end
  
    % if fieldmaporient is empty, then use interaction to get the values
    if isempty(fieldmaporient)
      figure; setfigurepos([100 100 800 400]);
      fieldmaporient = [0 0 0];
      while 1
        fun = @(x) flipdims(rotatematrix(x,1,2,fieldmaporient(1)),fieldmaporient(2:3));
        if isempty(inplanes)
          sltemp = round(linspace(1,size(epis{1},3),6)); sltemp = sltemp(2:5);
          subplot(1,2,1); imagesc(makeimagestack(epis{1}(:,:,sltemp))); colormap(gray); title('EPI');
        else
          sltemp = round(linspace(1,size(inplanes{1},3),6)); sltemp = sltemp(2:5);
          subplot(1,2,1); imagesc(makeimagestack(inplanes{1}(:,:,sltemp))); colormap(gray); title('Inplane');
        end
        sltemp = round(linspace(1,size(fieldmapbrain,3),6)); sltemp = sltemp(2:5);
        subplot(1,2,2); imagesc(makeimagestack(feval(fun,fieldmapbrain(:,:,sltemp)))); colormap(gray); title('Fieldmap brain');
        temp = input('We need to match the orientation of the fieldmap to the orientation of the other data.\nTo do this, we can apply CCW rotations and then flip the first and second dimensions as necessary.\nPlease specify a 3-element vector [A B C] where A is the number of CCW rotations to apply, B is whether to flip the first dimension, and C is whether to flip the second dimension (e.g. [1 1 0]).\nJust press RETURN if you are happy with how the current results are.\n--> ');
        if isempty(temp)
          fprintf('OK, the final value for <fieldmaporient> is %s.\n',mat2str(fieldmaporient));
          break;
        end
        fieldmaporient = temp;
      end
    end
  
    % rotate and flip based on fieldmaporient
    fun = @(x) flipdims(rotatematrix(x,1,2,fieldmaporient(1)),fieldmaporient(2:3));
    fieldmap = feval(fun,fieldmap);
    fieldmapbrain = feval(fun,fieldmapbrain);
  
    % record it
    fieldmaps{p} = fieldmap;
    fieldmapbrains{p} = fieldmapbrain;
    clear fieldmap fieldmapbrain;
    
  end
  
  % massage
  if ~iscell(fieldmapsizes)
    fieldmapsizes = {fieldmapsizes};
  end

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
  epis,episizes{1},epiinplanematrixsizes{1},cell2mat(epitr),episliceorder, ...
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
    save_nii(settr_nii(make_nii(int16(epis{p}),finalepisize),epitr{p}),sprintf(savefile,p));
  end
end
if ~isempty(savefileB)
  mkdirquiet(stripfile(savefileB));
  save_nii(make_nii(int16(validvol),finalepisize),savefileB);
end
if ~isempty(savefileC)
  mkdirquiet(stripfile(savefileC));
  save_nii(make_nii(int16(meanvol),finalepisize),savefileC);
end
fprintf('done (saving data).\n');

reportmemoryandtime;
