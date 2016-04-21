function [means, variances, weights, minLogL, posteriors, Si, expectation] = fitGMMwithSV(...,
    data, ...
    means, ...
    variances,...
    weights, ...
    theta, ...
    updateVariance, ...
    updateSpatialVariance, ...
    updateWeights, ...
    sameVariance, ...
    maxIters, ...
    convCriterion, ...
    trustWeights, ...
    trustVariance, ...
    numBins, ...
    verbosity, ...
    gmmFig, ...
    costFig)
% Description: Fits a Gaussian mixture model to data by means of Expectation Maximization
%
% Call: [means, variances, weights, minLogL, Si, dataEstimate] = fitGMM(data, means, variances, weights, useFixedMeans,
%           useSameVariance, maxMMIters, convCriterion, numBins, verbosity, figToPlotTo)
%
% Input:    data:               Data to fit the GMM to
%           prior:              Prior on tissue or structure classes
%           means:              Means of Gaussians in mixture
%           variances:          Variances of Gaussians in mixture
%           weights:            Weights of Gaussians in mixture
%           usedFixedMeans:     Means are updated according to the min and max value of the data (hence 2 parameters)
%           useSameVariance:    Equal variance is computed for all Gaussians
%           maxIters:           Maximum number of iterations allowed for fit
%           convCriterion:      Convergence threshold (relative change in cost)
%           numBins:            bins to plot histogram with
%           verbosity:          Verbosity level, higher value provides more output to stdout
%           figToPlotTo:        Figure handle to plot mixture model fitting to as iterations progress
%
% Output:   means:              Gaussian means that optimize likelihood of data
%           variances:          Gaussian variances that optimize likelihood of data
%           weights:            Gaussian weights that optimize likelihood of data
%           minLogL:            The negative loglikelihood at convergence
%           Si:                 The weights (or confidence) of each voxels in fit (posterior probability over variance)
%           dataEstimate:       The estimated data, given the optimized parameters
%
% Author: Christian Thode Larsen, christian.thode.larsen@gmail.com

% Copyright (c) 2015 Christian Thode Larsen, christian.thode.larsen@gmail.com
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

% optional parameter initialization
if nargin < 5
    theta = 1;
end
if nargin < 6
    updateVariance = 1;
end
if nargin < 7
updateSpatialVariance = 0;
end
if nargin < 8
    updateWeights = 1;
end
if nargin < 9
    sameVariance = 1;
end
if nargin < 10
    maxIters = 1000;
end
if nargin < 11
    convCriterion = 1e-5;
end
if nargin < 12
    trustWeights = 1e-3;
end
if nargin < 13
    trustVariance = 1e-3;
end
if nargin < 14
    numBins = 1024;
end
if nargin < 15
    verbosity = 'none';
end
if nargin < 16
    gmmFig = -1;
end
if nargin < 17
    costFig = -1;
end

% initialize
numChannels = size(means, 2);
numGaussians = length(means);
numVoxels = size(data,1);
posteriors = zeros(numVoxels, numGaussians);

% history of minLogLikelihoods - used to determine convergenge between iterations
costHist = 1/eps;
converged = false;

if ishandle(gmmFig) && ~sameVariance
    gaussianEvals = zeros(numBins, numGaussians);
    
    [histogram, binCenters] = hist(data(:,1), numBins);
    histogram = histogram / sum(histogram);
    binWidth = binCenters(2) - binCenters(1);
end

% Performing EM
for mmIter = 1:maxIters
    %plot the state of affairs
    if ishandle(gmmFig) && ~sameVariance
        figure(gmmFig)
        for k = 1:numGaussians
            gaussianEvals(:, k) =  (1 / sqrt(2*pi*variances(k, 1))) * ...
                exp(-0.5 * (binCenters - means(k, 1)).^2 / variances(k, 1)) * weights(k);
        end
        gaussianEvals = gaussianEvals * binWidth;
        
        clf
        grid
        bar(binCenters, histogram);
        hold on
        plot(binCenters, gaussianEvals', 'lineWidth', 1);
        plot(binCenters, sum(gaussianEvals, 2), 'r', 'lineWidth', 2);
        xlabel('Bins')
        xlim([1,7]);
        title('Log-intensity histogram');
        ylabel('Normalized Voxel Count');
        set(gca, 'fontsize', 6)
        axis tight
        hold off
        drawnow
    end
    
    for k = 1:numGaussians
        posteriors(:, k) = ((1 / sqrt(2*pi*variances(k,1))) * exp(-0.5 * ((data(:,1) - means(k,1)).^2) / (variances(k,1))));      
        for j = 2:numChannels
            posteriors(:, k) = posteriors(:, k).* ((1 / sqrt(2*pi*theta*variances(k,j))) * exp(-0.5 * ((data(:,j) - means(k,j)).^2) / (theta*variances(k,j))));
        end   
       posteriors(:,k) = posteriors(:,k) * weights(k);
    end
    
    normalizer = sum(posteriors, 2);
 
    for i = 1:numGaussians
        posteriors(:,i) = posteriors(:,i) ./ (normalizer + eps);
    end
    
    minLogL = -sum(log(normalizer + eps));
    costHist(end + 1) = minLogL;%#ok
    
    
    relChange = (costHist(end - 1) - costHist(end)) / numVoxels;
    if abs(relChange) < convCriterion || maxIters == 1
        if ~isempty(strfind(verbosity, 'debug'))
            fprintf('\tGGM fit converged in %i iterations\n', mmIter-1);
        end
        converged = 1;
        break
    end
    
    if ishandle(costFig)
        figure(costFig)
        plot(costHist(2:end));
        axis tight
        title('Cost');
        xlabel('Iteration');
        ylabel('Cost');
        drawnow
    end
    
    posteriorsSum = sum(posteriors);
    posteriorsSumSum = sum(posteriorsSum);
    
    if updateWeights
        % update weights - ' is necessary to maintain column-major
        weights = posteriorsSum' / posteriorsSumSum;
    end
    
     % don't update gaussians that contribute less than trustCriterion to the posterior ( which should == 1)
    trusted = find(weights >  trustWeights);
    if(~all(trusted))
        warning('Classes are collapsing');
    end
    
    % update means
    % made loop because of memory at full res
    newmeans = means;  
    for k = 1:numGaussians
        for j = 1:numChannels          
            newmeans(k,j) = (posteriors(:,k)' * data(:,j)) ./ (posteriorsSum(k) + eps);          
        end
    end
    means(trusted,:) = newmeans(trusted,:);
    
    if updateVariance
        % same variance for intensities!
        if sameVariance
            variance = 0;
            for k = 1:numGaussians
                variance = variance + (posteriors(:,k)' * (data(:,1) - means(k,1)).^2);
            end
            variance = variance / posteriorsSumSum;
            variances = repmat(variance, [numGaussians, 1]);
        % only update intensity variance, spatial variance is fixed
        else     
            newvariances = variances;
            for k = 1:numGaussians
                % update variances
                newvariances(k, 1) = (posteriors(:,k)' * (data(:,1) - means(k, 1)).^2) ./ (posteriorsSum(k) + eps);
                newvariances(k, 1) = max(newvariances(k, 1), trustVariance);

            end
            variances(trusted,1)= newvariances(trusted,1);
        end
    end
    if updateSpatialVariance
        if sameVariance
            for j = 2:numChannels
                variance = 0;
                for k = 1:numGaussians
                    variance = variance + (posteriors(:,k)' * (data(:,j) - means(k,j)).^2);
                end
                variance = variance / posteriorsSumSum;
                variances(:,j) = repmat(variance, [numGaussians, 1]);
            end
        else         
            newvariances = variances;
            for k = 1:numGaussians
                for j = 1:numChannels 
                    % update variances
                    newvariances(k,j) = (posteriors(:,k)' * (data(:,j) - means(k, j)).^2) ./ (posteriorsSum(k) + eps);
                    newvariances(k,j) = max(newvariances(k,j), trustVariance);
                end
            end
            variances(trusted,2:4)= newvariances(trusted,2:4);
        end
    end

end

if (~converged && ~isempty(strfind(verbosity, 'debug')))
    disp(['GMM stopped: Maximum number (', num2str(maxIters), ') of iterations reached'])
end

if(nargout > 6)
    Si = zeros(numVoxels, 1);
    expectation = zeros(numVoxels, 1);
    for k = 1:numGaussians
        Sik = posteriors(:,k) / variances(k,1);
        
        Si = Si + Sik;
        expectation = expectation + Sik * means(k,1);
    end
    expectation = expectation ./ Si;
end