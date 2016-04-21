function [mu,sigma] = compute_fwhm(input, num_bins)
% ------------------------------------------------------------------------
% [mu, sigma] = compute_fwhm(input, num_bins) 
% 
% Compute the Full Width Half Maximum FWHM of a voxel intensity tissue
% distribution. 
%  -input: 1D matrix containing the voxel intensities of a tissue
%    distribution
%  -num_bins: Number of bins used to generate the histogram.
%  
%  The function returns the \mu and \sigma as extracted from the generated
%  histogram
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
    
    [h_frequencies, h_centers] = hist(input, num_bins);    
    h_frequencies(1) = 0;
        
    [maxbin, pos_maxbin] = max(h_frequencies);
    
    % mu estimation: 
    mu = h_centers(pos_maxbin);
    
    % sigma estimation
    % lower band
    pos_hx1 = find(h_frequencies ./ maxbin > 0.5, 1,'first') -1;
    hx1 = h_frequencies(pos_hx1) / maxbin;
    hx1m1 = h_frequencies(pos_hx1 +1) / maxbin;
    
    mx1 = h_centers(pos_hx1);
    mx1m1 = h_centers(pos_hx1 +1);
    loband = mx1+(0.5-hx1)*(mx1m1-mx1)/(hx1m1-hx1);
    
    
    % higher band
    pos_hx2 = find((h_frequencies ./ maxbin) > 0.5, 1, 'last') +1;
    hx2 = h_frequencies(pos_hx2) / maxbin;
    hx2m1 = h_frequencies(pos_hx2 -1) / maxbin;
    
    mx2 = h_centers(pos_hx2);
    mx2m1 = h_centers(pos_hx2 -1);
    hiband = mx2m1 -(0.5-hx2)*(mx2m1-mx2)/(hx2m1-hx2);
    
    
    % sigma
    sigma=(hiband-loband)/(2*sqrt(2*log(2)));
   
 end