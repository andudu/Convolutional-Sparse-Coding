function im = stdimage(imnm, rflg, sbdr)
    
% stdimage -- Get a named standard image
%
% Usage:
%       im = stdimage(imnm, rflg, sbdr)
%
% Input:
%       imnm        String containing the image name
%       rflg        If true (default), read image and return data, otherwise
%                   return full path to image file
%       sbdr        Cell array of subdirectories to search (this
%                   argument should usually not be specified)
%
% Output:
%       im          Image data or image file path
%
%   
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2014-10-15
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


if nargin < 3.
  sbdr = {'Std', 'Kodak'};
end
if nargin < 2,
  rflg = 1;
end

% Determine base path to image data
p0 = which('sporco');
K = strfind(p0, '/');
p1 = p0(1:K(end)-1);
bp = [p1 '/Data'];

% Try to find specified file
ip = [];
if strcmp(imnm(end-3:end), '.png')
    imnm = imnm(1:end-4);
end

for k = 1:length(sbdr),
  ipt = [bp '/' sbdr{k} '/' imnm '.png' ];
  if exist(ipt,'file'),
    ip = ipt;
    break;
  end
end

if rflg && ~isempty(ip),
  im = imread(ip);
else
  im = ip;
end

return
