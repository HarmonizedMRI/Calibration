function [ims, b0, imsos1] = loadb04ge(pfile, deltaTE, readoutFile)
% function [ims, b0, imsos1] = loadb04ge(pfile, deltaTE, readoutFile)
%
% Reconstruct individual coil images and B0 map
%
% Inputs:
%  pfile         string
%  deltaTE       echo time offsets (ms), see b0.m
%  readoutFile   string. Default: 'readout.mod' inside 'b0.tar'
%
% Outputs:
%  b0            [nx ny nz] B0 field map calculated from the first two
%                entries in deltaTE (Hz)
%  ims           [nx ny nz ncoils] Individual coil images

if nargin < 3
    readoutFile = 'readout.mod';
    system('tar xf b0.tar readout.mod');
end

% Reconstruct individual coil images.
% 'recon3dft' assumes that readout.mod was created with toppe.utils.makegre()
for iecho = 1:length(deltaTE)
    ims(:,:,:,:,iecho) = toppe.utils.recon3dft(pfile, 'echo', iecho, ...
        'readoutFile', readoutFile, ...
        'flipFid', false);
end

% Calculate B0 map from first two echoes
if length(deltaTE) > 1
    phs = toppe.utils.phasecontrastmulticoil(ims(:,:,:,:,2), ims(:,:,:,:,1));
    b0 = phs/(2*pi)/(deltaTE(2)-deltaTE(1))*1e3;  % field map (Hz)
end

% root sum of squares coil-combined image from first echo
imsos1 = sqrt(sum(abs(ims(:,:,:,:,1)).^2, 4));

mask = imsos1 > 0.1*max(imsos1(:));
b0(~mask) = 0;

im(b0, [-100 100]); colormap default; colorbar;

