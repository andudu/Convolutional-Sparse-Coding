% This script sets up the Matlab path for the SPORCO library

% Add selected subdirectories of the parent directory of this
% script into the path
sporco_p0 = which('sporco');
sporco_K = strfind(sporco_p0, '/');
sporco_p1 = sporco_p0(1:sporco_K(end)-1);
sporco_sd = {'.', 'Demo', 'DictLearn', 'SparseCode', 'Util','Denoising',...
    'MultiDict','Nofilter','DictLReg','BM3D','BM3D/BM3D-SAPCA',...
    'BM3D/IDDBM3D','ompbox10','ksvdbox13','ZeroFill','SaltPepper',...
    'SuperRes','1dsignals','PatchComp'};
for sporco_k=1:length(sporco_sd),
  addpath([sporco_p1 '/' sporco_sd{sporco_k}]);
end
sporco_p2 = genpath([sporco_p1 '/Extrnl']);
addpath(sporco_p2);
sporco_path = sporco_p1;
clear sporco_p0 sporco_p1 sporco_p2 sporco_sd sporco_K sporco_k
