function outImage=intCorrect(image, mask)
% This function correct the intensities of the image via an anisotropic
% filter (by Perona and Malik) and N3 correction (by Thode).
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2016  Eloy Roura / Arnau Oliver / Xavier Llado
%   $Revision: 2.1$     $Date: 27/03/16$ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    im=load_untouch_nii(image);
    
    %Perona and Malik
       num_iter = 1;
       delta_t = 3/44;
       kappa = 50;
       option = 1;
       voxel_spacing = ones(3,1);
       im_denoised = anisodiff3D(im.img, num_iter, delta_t, kappa, option, voxel_spacing);
       
    
    %N3 correction    
    im_denoised(~mask)=0;    
    mask(im_denoised==0)=0;
    
    [outImage, biasField, maskOut] = iic(double(im_denoised), ...
        'smoothingDistance', 30, ...
        'stepSize', [im.hdr.dime.pixdim(2), im.hdr.dime.pixdim(3), im.hdr.dime.pixdim(4)], ...
        'smoothingRegularization', 1e-4, ...
        'backgroundThresholdType','otzu');
        %'mask', double(mask));
   
end