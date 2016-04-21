function [lesion_mask, thresholded_mask] = lesionSegmentationTool(input_image, tissueProb, alpha, omega_t, omega_n, min_size, ventricles)
% ------------------------------------------------------------------------
% [lesion_mask] = lesionSegmentationTool(input_image, tissueProb, alpha, omega_t, omega_n, min_size) 
% 
% WM lesion segmentation of MS lesions on FLAIR images. 
%
%
%  -input_image  --> Pre-processed FLAIR (skull-stripped + intensity corrected) 
%  -tissueProb   --> labeled T1-w tissue segmentation (1 = CSF, 2 = GM, 3 = WM)
%  -alpha:       --> Parameter for controlling the weight of the sigma computed from 
%                    [mu, sigma] = FWHM function.
%  - omega_t     --> Percentage of either GM / WM for candidate lesion regions
%  - omega_n     --> Percentage of WM voxels for 6-connectivity 3D boundary
%                    voxels of candidate_regions.
%   -min_size    --> minimum area of candidate lesion regions.
%   -ventricles  --> ventricles mask to remove its attached elongated lesions.
%
%  - lesion_mask  --> Binary mask of the final candidate lesion regions.
%  - thresholded_mask  --> Binary mask of the initial candidate lesion regions.
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
%   Copyright (C) 2016  Eloy Roura / Sergi Valverde / Arnau Oliver / Xavier Llado
%   $Revision: 2.1$     $Date: 27/03/16$ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % vars:
    lesion_mask = zeros(size(input_image)); % output lesionmask
    gm_mask = tissueProb == 2;              % GM TISSUE
    wm_mask = tissueProb == 3;              % WM TISSUE

    
    % ********************************************************************
    % 1. Threshold evaluation
    %
    % - by default: 512 bins are computed. The final threshold is computed
    %   following T = mu + alpha * sigma.
    % - So far, the FWHM methods are different between implementations :(
    %
    % ********************************************************************
    
    gm_flair = input_image(gm_mask == 1);
    [mu, sigma] = compute_fwhm(gm_flair,512);
    T = mu + (alpha * sigma);
    thresholded_mask = input_image >= T;
    fprintf('\t\tThreshold: %s + (%s * %s)\n', num2str(T),num2str(alpha),num2str(sigma));
    
    %disp(['Threshold: ', num2str(T), ' (',num2str(mu),' + ',num2str(alpha),' * ',num2str(sigma),')']);
    
    
    % ********************************************************************
    % 2. Filtering step
    %
    % - Connected components are computed using 6-neighbor connectivity in
    %   3D. (followin c++ implementation)
    % ********************************************************************
    
    % Find connected components in 3D
    CC = bwconncomp(thresholded_mask,6);
    candidate_regions = labelmatrix(CC);
    candidate_labels = unique(nonzeros(candidate_regions));
    
    fprintf('\t\tNumber of lesions before refinement: %s\n', num2str(CC.NumObjects));
    %disp(['Number of lesions before refinement: ', num2str(CC.NumObjects)]);
    
    % RULE 1). 3d lesion size has to be higher than min_size parameter
    labels_filter_1 = cellfun(@(x) numel(x) < min_size, CC.PixelIdxList);
    candidate_labels(labels_filter_1) = 0;
    
    
    % RULE 2). omega_t percentage of lesion voxels either GM and WM 
    labels_filter_2 = cellfun(@(x) ((sum(gm_mask(x)) + sum(wm_mask(x))) / numel(x)) < omega_t, CC.PixelIdxList);
    candidate_labels(labels_filter_2) = 0;  
  
    
    % RULE 3) Omega_n Percentage of neighbors as WM
    % so far, have to visit each lesion candidate with a for loop :(
    candidate_labels = nonzeros(candidate_labels);
    num_lesions = 0;
    for l=candidate_labels'
            % current label 
            region_voxels = find(candidate_regions == l);
            
            % default config for neighbors: 24 connectivity in 3D (radius 1)
            % Neighbor voxels for each region are found by substracting the
            % total neighbors of the candidate region and the same
            % candidate region. Not shure if efficient :/
            boundary_region_voxels = compute_neighborhoods(region_voxels, size(candidate_regions), 1,3);
            candidate_regions(boundary_region_voxels) = l;
            candidate_regions(region_voxels) = 0;
            boundary_voxels = find(candidate_regions == l);
            
            %boundary_voxels = boundary_region_voxels{1};
            perc_wm_elem = sum(wm_mask(boundary_voxels)) / numel(boundary_voxels);
            if perc_wm_elem >= omega_n
                lesion_mask(region_voxels) = 1;
                num_lesions = num_lesions +1;
            end
    end
    
    
    %New RULE. Aviod elongated lesions attached to the ventricles.
    if(~strcmp(ventricles,'none'))
     CC = bwconncomp(lesion_mask, 18);
     candidate_regions = labelmatrix(CC);
     candidate_labels = unique(nonzeros(candidate_regions));
     
     num_lesions = max(candidate_labels);
     for l=candidate_labels'
        region_voxels = find(candidate_regions == l);
        
        %If region_voxels intersects the ventricles
        if find(ventricles(region_voxels)==1)>1
            [x,y,z]=ind2sub(size(input_image),region_voxels);
            bb_voxels2D=(max(x)-min(x))*(max(y)-min(y))*(max(z)-min(z));
 
            if numel(region_voxels)>100 && ((max(x)-min(x))/(max(y)-min(y))<0.4 || (max(y)-min(y))/(max(x)-min(x))<0.4 || numel(region_voxels)<(bb_voxels2D*0.5)) 
                lesion_mask(region_voxels) = 0;
                num_lesions = num_lesions -1;
            end
        end       
     end
    end
    
    fprintf('\t\tFinal number of lesions: %s\n\n', num2str(num_lesions));
end