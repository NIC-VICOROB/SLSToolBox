% *****************************************************
%
% SPM12 batch code
%
% function call: matlab -nodisplay -nosplash -nodesktop -nojvm -r "deformBack('rT1.nii', 'y_rT1.nii'); quit"
% *****************************************************

function deformBack(path, defField, structure)
%
% This is an SPM batch code to deform one image (structure) using a 
% deformation field (defField) previously computed. The output will be
% stored in path.
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

    clc;
    spm_jobman('initcfg');

    matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def = {[path defField]};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[spm('dir'),'/toolbox/SLS/structures/', structure]};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    
    spm('defaults','pet');    
    spm_jobman('run',matlabbatch);    
    spm_clf;
    clear matlabbatch;
    
    %system(['mv ' [spm('dir'),'/toolbox/SLS/structures/w' structure ' '] path '/w' structure ' ']);
    movefile(['w' structure], [spm('dir'),'/toolbox/SLS/structures/w' structure]);

end