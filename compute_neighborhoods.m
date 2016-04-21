function [neighbor_positions] = compute_neighborhoods(target_positions, original_size, n, neighbors_dimension)
% ------------------------------------------------------------------------
% [neighbor_voxels] = compute_neighborhoods(target_positions, original_size, n, neighbors_dimension) 
% 
% Compute the positions of the {n x n x n} neighbors voxel given a  list of target
% positions. 
%
%  -target_positions: vector containing the positions of the 3D matrix
%                     expressed as linear indices
%  -original_size: Size of the original 3D matrix
%  -n: radius of the number of neighbors (total neigh: (2*n +1)^2 in 2D
%      or (2*n+1)^3 for 3D matrices.
%  -neighbors_dimension: Compute either the neighbor positions in 2D or 
%                        3D.   
%
%  neighbor_voxels: returns a matrix with size [num voxels, n^neighbor_dimension] 
%                   with the positions ocompute_neighborhoods3.mf all the adjacent neighbors.
%
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2016  Sergi Valverde / Eloy Roura / Arnau Oliver / Xavier Llado
%   $Revision: 2.1$     $Date: 27/03/16$ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if neighbors_dimension ~= 2 && neighbors_dimension ~= 3
        error([num2str(neighbors_dimension), ' appears not a valid neighbor dimensionality argument']);
    end
    
    diameter = (2*n+1);
    
    original_r = original_size(1);
    original_c = original_size(2);
    original_s = original_size(3);
    [rv,cv,sv] = ind2sub(original_size, target_positions); 

    % compute based on expanded matrix. We generate all the
    % neighbor positions by expanding them (exclude the center voxel)
    expansor = repmat(-n:n, numel(rv),1);
    N_r = repmat(rv,1,diameter) + expansor;
    N_c = repmat(cv,1,diameter) + expansor;

    % check dimensionality 
    if neighbors_dimension == 2
        N_s = repmat(sv,1,diameter); 
        % rows, cols and slices are concatenated in the column dimension to
        % generate all possible neighbor combinations for each voxel
        expanded_rows = repmat(N_r,1, diameter);  
        expanded_cols = kron(N_c, ones(1,diameter)); 
        expanded_slices = kron(N_s, ones(1,diameter));
        neighbor_positions = sub2ind(original_size, max(1, min(original_r, expanded_rows)), ...
                                                    max(1, min(original_c, expanded_cols)), ...
                                                    max(1, min(original_s, expanded_slices)));

    else
        N_s = repmat(sv,1,diameter) + expansor; 
        expanded_rows = repmat(N_r,1, diameter^2);  
        expanded_cols = repmat(kron(N_c, ones(1,diameter)),1,diameter); 
        expanded_slices = kron(N_s, ones(1,diameter^2));
         neighbor_positions = sub2ind(original_size, max(1, min(original_r, expanded_rows)), ...
                                                    max(1, min(original_c, expanded_cols)), ...
                                                    max(1, min(original_s, expanded_slices)));
    end
end 
