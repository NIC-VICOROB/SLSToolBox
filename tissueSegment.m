function [segmTissue, brainMask]=tissueSegment(input_folder, input_volume)
% This function is based on the SPM12 segment. The output of this function
% consists of the 3 main tissue maps saved in three different files in the  
% input_folder.
%
% Besides, after the spm tissue segmentation, we perform a maximum likelihood
% thresholded at 0.5 to obtain the 3 tissue masks and the brain mask:
%
%       BrainMask = P(WM|GM|CSF)>0.5;
%       TissueMask = Max(WM,GM,CSF);
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

    [~, spm_version, ~] = spm_fileparts(spm('dir'));
    if strcmp(spm_version,'spm12')
        TMP_path = [spm('dir'),'/tpm/TPM.nii'];
        p_tools = 'spatial';
        p_preproc = 'preproc';    
    else 
        TMP_path = [spm('dir'),'/toolbox/Seg/TPM.nii'];
        p_tools = 'tools';
        p_preproc = 'preproc8';
    end

	spm_jobman('initcfg');

	matlabbatch{1}.spm.(p_tools).(p_preproc).channel.vols = {[input_folder '/' input_volume]};
	matlabbatch{1}.spm.(p_tools).(p_preproc).channel.biasreg = 0.0001;
	matlabbatch{1}.spm.(p_tools).(p_preproc).channel.biasfwhm = 60;
	matlabbatch{1}.spm.(p_tools).(p_preproc).channel.write = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(1).tpm = {[TMP_path,',1']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(1).ngaus = 2;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(1).native = [1 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(1).warped = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(2).tpm = {[TMP_path,',2']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(2).ngaus = 2;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(2).native = [1 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(2).warped = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(3).tpm = {[TMP_path,',3']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(3).ngaus = 2;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(3).native = [1 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(3).warped = [0 0];	
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(4).tpm = {[TMP_path,',4']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(4).ngaus = 3;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(4).native = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(4).warped = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(5).tpm = {[TMP_path,',5']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(5).ngaus = 4;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(5).native = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(5).warped = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(6).tpm = {[TMP_path,',6']};
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(6).ngaus = 2;
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(6).native = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).tissue(6).warped = [0 0];
	matlabbatch{1}.spm.(p_tools).(p_preproc).warp.reg = 4;
	matlabbatch{1}.spm.(p_tools).(p_preproc).warp.affreg = 'mni';
	matlabbatch{1}.spm.(p_tools).(p_preproc).warp.samp = 3;
	matlabbatch{1}.spm.(p_tools).(p_preproc).warp.write = [1 0];

	spm_jobman('run',matlabbatch);
    clear matlabbatch;
    
    %Masking the probability maps
    c1 = load_untouch_nii([input_folder, '/c1', input_volume]);
    c2 = load_untouch_nii([input_folder, '/c2', input_volume]);
    c3 = load_untouch_nii([input_folder, '/c3', input_volume]);
    
    c1.hdr.dime.scl_slope=1;
    c2.hdr.dime.scl_slope=1;
    c3.hdr.dime.scl_slope=1;

    brainMask = ones(size(c1.img));
    brainMask(c1.img<0.5 & c2.img<0.5 & c3.img<0.5)=0;
    seg = zeros(size(c1.img));
    [~,seg] = max([c1.img(brainMask==1) c2.img(brainMask==1) c3.img(brainMask==1)],[],2);
    
    res=zeros(size(c1.img));
    res(brainMask==1)=seg;
    segmTissue=zeros(size(c1.img));
    segmTissue(res==3)=1;
    segmTissue(res==1)=2;
    segmTissue(res==2)=3;

    disp(['Tissue segmentation on volume ', input_volume, ' done!']);
end