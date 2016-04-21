function lesionSeg_run(job)
% SLS segments WM lesions wich are outliers of the GM distribution.
%
% SLS(inputT1_scan, inputFLAIR_scan, brain_mask, options, params) saves 
% a NIFTI  tissue map segmentation of the T1 scan and lesion mask 
% segmentation of the FLAIR scan. Output directories are the same as the 
% input images.
% 
%   tissue maps = pathT1/spmSegmTissue.nii
%   lesion mask = pathFLAIR/segmLesions.nii
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

%% Initial variables
	
    [pthT1, namT1, extT1] = spm_fileparts(job.data_T1{1});
    T1 = fullfile(pthT1, [namT1, extT1]);
    T1image=load_untouch_nii(T1);
    
    [pthFLAIR, namFLAIR, extFLAIR] = spm_fileparts(job.data_FLAIR{1});
    FLAIR = fullfile(pthFLAIR, [namFLAIR, extFLAIR]);
    FLAIRimage=load_untouch_nii(FLAIR);
    
    im2write=FLAIRimage;


    %% Tissue segmentation by SMP8/12
    
    if(max(T1image.img(:))>3)
        
        coreg(FLAIR, pthT1, [namT1, extT1]);
        [segmTissues, brainMask]=tissueSegment(pthT1, ['r' namT1]);
        im2write.img=segmTissues;
        save_untouch_nii(im2write, [pthT1, '/spmSegmTissues.nii']);
        
        if(job.exclusion.ventricles)
            deformBack(pthT1, ['/iy_r' namT1 '.nii'], 'ventricles.nii');
        end
            
        %Clean up tissue segmentation files
        if (job.outputs.cleanup)
            delete([pthT1, '/c*' namT1, extT1]);
            delete([pthT1, '/*seg8.mat']);
            delete([pthT1, '/*seg8.txt']);
        end
        
        
        FLAIRimage.img=intCorrect(FLAIR, brainMask);
        save_untouch_nii(FLAIRimage, [pthT1,'/FLAIR.corrected.nii']);
        
    else
        segmTissues = T1image.img;
        FLAIRimage = load_untouch_nii(FLAIR);
    end
    
    
    
    %% Lesion segmentation method by binary file
    
    %%1st iteration
    fprintf('\n\t1st iteration:\n');
    
	if(job.params.stIter.omegaT < 0 || job.params.stIter.omegaT > 1)
	    error('Accepted values for omegaT are [0...1].')
	end

	if(job.params.stIter.omegaN < 0 || job.params.stIter.omegaN > 1)
	    error('Accepted values for omegaN are [0...1].')
	end

	if(job.params.stIter.lesionSize < 0 )
	    error('Lesion Size must be a positive integer.')
	end
	
	alpha = job.params.stIter.alpha;
	omegaT = job.params.stIter.omegaT;
	omegaN = job.params.stIter.omegaN;
	lesionSize = job.params.stIter.lesionSize;
    
    
    if(job.exclusion.ventricles)
        ventricles=load_untouch_nii([pthT1 '/wventricles.nii']);
        window=strel('disk',3);
        ventricles.img(ventricles.img~=1)=0;
        ventricles=imdilate(ventricles.img,window);
    else
        ventricles='none';
    end
    
    [lesion_mask1, thr1] = lesionSegmentationTool(double(FLAIRimage.img), segmTissues, alpha, omegaT, omegaN, lesionSize, ventricles);
    
	FLAIRimage.img(lesion_mask1==1)=0;
    segmTissues(lesion_mask1==1)=0;
    
    %%2nd iteration
    fprintf('\t2nd iteration:\n');
        
    if(job.params.ndIter.omegaT2 < 0 || job.params.ndIter.omegaT2 > 1)
        error('Accepted values for omegaT2 are [0...1].')
    end

    if(job.params.ndIter.omegaN2 < 0 || job.params.ndIter.omegaN2 > 1)
        error('Accepted values for omegaN2 are [0...1].')
    end

    if(job.params.ndIter.lesionSize2 < 0 )
        error('Lesion Size must be a positive integer.')
    end
	
    alpha2 = job.params.ndIter.alpha2;
    omegaT2 = job.params.ndIter.omegaT2;
    omegaN2 = job.params.ndIter.omegaN2;
    lesionSize2 = job.params.ndIter.lesionSize2;
        
    [lesion_mask2, thr2] = lesionSegmentationTool(double(FLAIRimage.img), segmTissues, alpha2, omegaT2, omegaN2, lesionSize2, ventricles);

    lesion_mask=lesion_mask1;
    lesion_mask(lesion_mask2==1)=1;
    
    T1image.img=lesion_mask;
    
    save_untouch_nii(T1image, [pthFLAIR, '/segmLesions.nii']);
	
    
    fprintf('Lesion segmentation on volume %s done!\n\n', FLAIR);
    
    %% Clean up temporary files
   
    if (job.outputs.lesionMaskst)
        T1image.img=lesion_mask1;
        save_untouch_nii(T1image, [pthFLAIR, '/segmLesions.iter1.nii']);
    end
    
    if (job.outputs.lesionMasknd)
        T1image.img=lesion_mask2;
        save_untouch_nii(T1image, [pthFLAIR, '/segmLesions.iter2.nii']);
    end
    
     if (job.outputs.thrst)         
         T1image.img=thr1;
         save_untouch_nii(T1image, [pthFLAIR, '/thr.iter1.nii']);
     end
     
     if (job.outputs.thrnd)
         T1image.img=thr2;
         save_untouch_nii(T1image, [pthFLAIR, '/thr.iter2.nii']);
     end

     