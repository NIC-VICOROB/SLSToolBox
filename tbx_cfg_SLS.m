function SLS = tbx_cfg_SLS
%     SLS segments WM lesions wich are outliers of the GM distribution.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2016  Eloy Roura / Arnau Oliver / Xavier Llado
%   $Revision: 2.1$     $Date: 27/03/16$ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------------------------------------------------------
% Input Images for lesion and tissue segmentation
% ---------------------------------------------------------------------
data_T1 = cfg_files;
data_T1.tag  = 'data_T1';
data_T1.name = 'T1 volume';
data_T1.help = {'Select raw T1 image.'
                'We recomend to perform a skull stripping process, denoising and a bias correction beforehand, as expleined in the help menu.'
                'Image orientation into the ICBM/MNI space is also required.'
                ''
                'Instead of the T1w raw image tissue segmentation image can be selected.'
                ''
                '*ATTENTION! Tissue labels must follow the notation below:'
                '	Cerebrospinal fluid -> 1'
                '	Grey matter -> 2'
                '	White matter -> 3'
};
data_T1.filter = 'image';
data_T1.ufilter = '.nii';
data_T1.num     = [1 Inf];

data_FLAIR = cfg_files;
data_FLAIR.tag  = 'data_FLAIR';
data_FLAIR.name = 'FLAIR volume';
data_FLAIR.help = {'Select raw FLAIR image. We recomend to perform a skull stripping process, denoising and a bias correction beforehand, as expleined in the help menu.' 
                    'Furthermore, T1w image and FLAIR must be aligned and ICBM/MNI oriented.'};
data_FLAIR.filter = 'image';
data_FLAIR.ufilter = '.nii';
data_FLAIR.num     = [1 Inf];


% ---------------------------------------------------------------------
% alpha
% ---------------------------------------------------------------------
alpha         = cfg_entry;
alpha.tag     = 'alpha';
alpha.name    = 'Deviations over the GM mean (alpha)';
alpha.help    = {'Number of deviations over GM mean defining the lesions. The threshold is computed using the following formula:'
                '	Thr = mean + alpha*std'
                'where the mean and standard deviation is computed over the GM histogram of the FLAIR image.'
                'Recomendations:'                    
                    '	For 1 iteration: 2.5'
                    '	For 2 iteration: 3'
};
alpha.num     = [1 1];
alpha.val     = {2.5};

% ---------------------------------------------------------------------
% OmegaT
% ---------------------------------------------------------------------
omegaT         = cfg_entry;
omegaT.tag     = 'omegaT';
omegaT.name    = 'WM/GM/GM+CSF tissue % of the lesion (Lambda_ts)';
omegaT.help    = {'This parameter is used to define the percentage of voxels inside the lesion segmented as WM, GM or PV (GM+CSF) is required (??ts). Lesions with values below this threshold are discarded. Values between [0-1] are allowed where 0 is less restrictive than 1.'};
omegaT.num     = [1 1];
omegaT.val     = {0.6};

% ---------------------------------------------------------------------
% OmegaN
% ---------------------------------------------------------------------
omegaN         = cfg_entry;
omegaN.tag     = 'omegaN';
omegaN.name    = 'WM tissue % of the lesion neighborhood (lambda_nb)';
omegaN.help    = {'This parameter is used to define the percentage of lesion surrounding voxels segmented as WM voxels is required (??nb). Lesions with values below this threshold are discarded. Values between [0-1] are allowed where 0 is less restrictive than 1.'
                    'Recomendations:'                    
                    '	For 1 iteration: 0.4-0.6'
                    '	For 2 iterations: 1st = 0.6-0.7'
};
omegaN.num     = [1 1];
omegaN.val     = {0.55};

% ---------------------------------------------------------------------
% lesionSize
% ---------------------------------------------------------------------
lesionSize         = cfg_entry;
lesionSize.tag     = 'lesionSize';
lesionSize.name    = 'Lesion Size';
lesionSize.help    = {'Minimum size accepted as an MS lesion.'
			'Recomendations (30mm3):'
			'	1x1x3mm = 10 voxels (30mm3)'
			'	1x1x1mm = 30 voxels (30mm3)'
}';
lesionSize.num 	   = [1 1];
lesionSize.val     = {30};

% ---------------------------------------------------------------------
% alpha2
% ---------------------------------------------------------------------
alpha2         = cfg_entry;
alpha2.tag     = 'alpha2';
alpha2.name    = 'Deviations over the GM mean (alpha)';
alpha2.help    = {'Number of deviations over GM mean defining the lesions. The threshold is computed using the following formula:'
                '	Thr = mean + alpha*std'
                'where the mean and standard deviation is computed over the GM histogram of the FLAIR image.'
}';
alpha2.num     = [1 1];
alpha2.val     = {2};

% ---------------------------------------------------------------------
% OmegaT2
% ---------------------------------------------------------------------
omegaT2         = cfg_entry;
omegaT2.tag     = 'omegaT2';
omegaT2.name    = 'WM/GM/GM+CSF tissue % of the lesion (lambda_ts)';
omegaT2.help    = {'This parameter is used to define the percentage of voxels inside the lesion segmented as WM, GM or PV (GM+CSF) is required (??ts). Lesions with values below this threshold are discarded. Values between [0-1] are allowed where 0 is less restrictive than 1.'};
omegaT2.num     = [1 1];
omegaT2.val     = {0.75};

% ---------------------------------------------------------------------
% OmegaN2
% ---------------------------------------------------------------------
omegaN2         = cfg_entry;
omegaN2.tag     = 'omegaN2';
omegaN2.name    = 'WM tissue % of the lesion neighborhood (lambda_nb)';
omegaN2.help    = {'This parameter is used to define the percentage of lesion surrounding voxels segmented as WM voxels is required (??nb). Lesions with values below this threshold are discarded. Values between [0-1] are allowed where 0 is less restrictive than 1.'};
omegaN2.num     = [1 1];
omegaN2.val     = {0.7};

% ---------------------------------------------------------------------
% lesionSize2
% ---------------------------------------------------------------------
lesionSize2         = cfg_entry;
lesionSize2.tag     = 'lesionSize2';
lesionSize2.name    = 'Lesion Size';
lesionSize2.help    = {'Minimum size accepted as an MS lesion.'
			'Recomendations (20-30mm3):'
			'	1x1x3mm = 10 voxels (30mm3)'
			'	1x1x1mm = 30 voxels (20mm3)'
}';
lesionSize2.num 	   = [1 1];
lesionSize2.val     = {20};


% ---------------------------------------------------------------------
% Exclusion regions
% ---------------------------------------------------------------------
ventricles         = cfg_menu;
ventricles.tag     = 'ventricles';
ventricles.name    = 'Ventricles';
ventricles.help    = {'To Exclude candidates next to the ventricles.'}';
ventricles.labels  = {
               'yes'
               'none'
}';
ventricles.values  = {1 0};
ventricles.val     = {0};




% ---------------------------------------------------------------------
% Lesion mask first iter
% ---------------------------------------------------------------------
lesionMaskst         = cfg_menu;
lesionMaskst.tag     = 'lesionMaskst';
lesionMaskst.name    = 'Lesion mask 1st iteration';
lesionMaskst.help    = {'Save the lesion mask of the first iteration.'}';
lesionMaskst.labels  = {
               'yes'
               'none'
}';
lesionMaskst.values  = {1 0};
lesionMaskst.val     = {0};

% ---------------------------------------------------------------------
% Lesion mask second iter
% ---------------------------------------------------------------------
lesionMasknd         = cfg_menu;
lesionMasknd.tag     = 'lesionMasknd';
lesionMasknd.name    = 'Lesion mask of the 2nd iteration';
lesionMasknd.help    = {'Save the lesion mask of the second iteration, i.e. over the candidates of the second iteration.'};
lesionMasknd.labels  = {
               'yes'
               'none'
}';
lesionMasknd.values  = {1 0};
lesionMasknd.val     = {0};

% ---------------------------------------------------------------------
% Candidates 1st iter
% ---------------------------------------------------------------------
thrst         = cfg_menu;
thrst.tag     = 'thrst';
thrst.name    = 'Candidates 1st iteration';
thrst.help    = {'Save the candidate lesions of the first iteration'}';
thrst.labels  = {
               'yes'
               'none'
}';
thrst.values  = {1 0};
thrst.val     = {0};

% ---------------------------------------------------------------------
% Candidates 2nd iter
% ---------------------------------------------------------------------
thrnd         = cfg_menu;
thrnd.tag     = 'thrnd';
thrnd.name    = 'Candidates 2nd iteration';
thrnd.help    = {'Save the candidate lesions of the second iteration'}';
thrnd.labels  = {
               'yes'
               'none'
}';
thrnd.values  = {1 0};
thrnd.val     = {0};

% ---------------------------------------------------------------------
% Temporary files
% ---------------------------------------------------------------------
cleanup         = cfg_menu;
cleanup.tag     = 'cleanup';
cleanup.name    = 'Temporary files';
cleanup.help    = {'Save the temporary files'}';
cleanup.labels  = {
               'yes'
               'none'
}';
cleanup.values  = {1 0};
cleanup.val     = {0};






% ---------------------------------------------------------------------
% First iteration
% ---------------------------------------------------------------------
stIter         = cfg_branch;
stIter.tag     = 'stIter';
stIter.name    = 'First iteration';
stIter.val     = {alpha omegaT omegaN lesionSize};
stIter.help    = {'Parameters for the first iteration'};


% ---------------------------------------------------------------------
% Second iteration
% ---------------------------------------------------------------------
ndIter         = cfg_branch;
ndIter.tag     = 'ndIter';
ndIter.name    = 'Second iteration';
ndIter.val     = {alpha2 omegaT2 omegaN2 lesionSize2};
ndIter.help    = {'Parameters for the second iteration'};


% ---------------------------------------------------------------------
% Parameters
% ---------------------------------------------------------------------
params         = cfg_branch;
params.tag     = 'params';
params.name    = 'Parameters';
params.val     = {stIter ndIter};
params.help    = {'Parameters for lesion segmentation'};


% ---------------------------------------------------------------------
% Writting options
% ---------------------------------------------------------------------
outputs         = cfg_branch;
outputs.tag     = 'outputs';
outputs.name    = 'Output options';
%outputs.val     = {lesionMask lesionMaskst lesionMasknd thrst thrnd cleanup};
outputs.val     = {lesionMaskst lesionMasknd thrst thrnd cleanup};
outputs.help    = {'Files to write at the end of the process'};


% ---------------------------------------------------------------------
% Exclusion options
% ---------------------------------------------------------------------
exclusion         = cfg_branch;
exclusion.tag     = 'exclusion';
exclusion.name    = 'Exclusion options';
exclusion.val     = {ventricles};
exclusion.help    = {'Regions to exclude when removing potential candidates'};



% ---------------------------------------------------------------------
% Lesion Segmentation
% ---------------------------------------------------------------------
lesionSegment         = cfg_exbranch;
lesionSegment.tag     = 'lesionSegment';
lesionSegment.name    = 'Lesion segmentation tool';
lesionSegment.val     = {data_T1 data_FLAIR params exclusion outputs};
lesionSegment.help    = {'This toolbox is able to segment WM MS lesions over FLAIR images.'
    'The method implementd by this tool is described in Roura et al. (2015), which follows the steps below:'
    ''
    'Inputs:  - T1-w image or its tissue segmentation'
    '         - Flair image.'
    ''
    '1. Pre-processing'
        '  - Skull-stripping by spm tissue segmentation.'
        '  - Image denoising by the N-dimensional version of the classic Perona and Malik (1990).'
        '  - Bias correction by N3 Sled et al. (1998) (matlab implementation by Thode et al. (2014)).'
        '  - Intra-subject co-registration by SPM affine registration and ICBM/MNI orientation.'
    ''
    '2. Process'
        'Iter 1: bigger and higher intensity lesion segmentation -> ++lesion size, --restrictions.'
        'Iter 2: smaller and lower intensity lesion segmentation -> --lesion size, ++restrictions.'
  };
lesionSegment.prog    = @lesionSeg_run;


SLS  = cfg_choice;
SLS.name = 'SLS';
SLS.tag  = 'SLS';
SLS.values = {lesionSegment};
