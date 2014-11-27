% make directory for corrected DICOMs
mkdirquiet(newdatadir);

% loop over DICOM directories
for p=1:length(epifilenames)
  fprintf('processing DICOM directory %s\n',epifilenames{p});

  % make directory
  mkdirquiet(fullfile(newdatadir,stripfile(epifilenames{p},1)));

  % load corrected data
  fprintf('loading data from %s\n',sprintf(savefile,p));
  a1 = load_untouch_nii(sprintf(savefile,p));

  % figure out all the DICOM files
  b1 = matchfiles([epifilenames{p} '/*.dcm']);
  assert(size(a1.img,4)==length(b1));  % sanity check

  % process each file
  for q=1:length(b1)
    statusdots(q,length(b1));
  
    % load meta data
    metadata = dicominfo(b1{q});
    
    % make a mosaic
    nn = ceil(sqrt(size(a1.img,3)));
    extrapad = repmat({zeros(sizefull(a1.img,2),dataformat)},[1 nn*nn-size(a1.img,3)]);
    temp = cell2mat(reshape([splitmatrix(cast(a1.img(:,:,:,q),dataformat),3) extrapad],[nn nn])');
    
    % write corrected data to new file
    dicomwrite(temp, ...
               fullfile(newdatadir,stripfile(epifilenames{p},1),stripfile(b1{q},1)), ...
               metadata,'CreateMode','copy');
  
  end
  fprintf('\n');

end
