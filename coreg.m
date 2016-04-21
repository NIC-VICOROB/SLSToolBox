% *****************************************************
%
% SPM8 batch code
%
%function call: matlab -nodisplay -nosplash -nodesktop -nojvm -r "spm8_reg('T1.nii', 'FLAIR.nii'); quit"
% *****************************************************

function output=coreg(target, sourcePath, sourceImage)
%
% This is an SPM batch code to coregister two input images. The sourceImage
% is aligned to the targetImage and stored in sourcePath preceded by an the
% char 'r'.
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

    %clc;
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[target, ',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[[sourcePath '/' sourceImage], ',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    spm_jobman('run',matlabbatch);    
    clear matlabbatch;

    output=[sourcePath '/r' sourceImage];
end