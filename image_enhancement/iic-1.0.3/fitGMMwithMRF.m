function [means, variances, weights, gammas, minLogL, Si, expectation] = fitGMMwithMRF(...,
    data, ...
    means, ...
    variances,...
    weights, ...
    gammas, ...
    equidistantMeans, ...
    updateMeans, ...
    sameVariance, ...
    updateVariance, ...
    betaMRF, ...
    updatePatternMRF, ...
    maxIters, ...
    convCriterion, ...
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
%           means:              Means of Gaussians in mixture
%           variances:          Variances of Gaussians in mixture
%           weights:            Weights of Gaussians in mixture
%           usedFixedMeans:     Means are updated according to the min and max value of the data (hence 2 parameters)
%           useSameVariance:    Equal variance is computed for all Gaussians
%           betaMRF:            neighborhood regularization, such that a given voxel want to belong to the same voxel as its neighbors
%           updatePatternMRF:   updatePattern that can be precomputed for conditionally independent voxels
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
if nargin < 6
    equidistantMeans = 0;
end
if nargin < 7
    updateMeans = 1;
end
if nargin < 8
    sameVariance = 0;
end
if nargin < 9
    updateVariance = 1;
end
if nargin < 10
    betaMRF = 0;
end
if nargin < 11
    updatePatternMRF = [];
end
if nargin < 12
    maxIters = 100;
end
if nargin < 13
    convCriterion = 1e-5;
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
numGaussians = length(means);
numVoxels = numel(data);

squaredDists = zeros(numVoxels, numGaussians);
likelihoods = zeros(numVoxels, numGaussians);

alphas = 1 - (0:numGaussians-1)' / (numGaussians-1);

% history of minLogLikelihoods - used to determine convergence between iterations
costHist = 1/eps;
converged = false;

if ishandle(gmmFig)
    gaussianEvals = zeros(numBins, numGaussians);
    [histogram, binCenters] = hist(data, numBins);
    histogram = histogram / sum(histogram(:));
    binWidth = binCenters(2) - binCenters(1);
end

% Performing EM
for mmIter = 1:maxIters  
    if ishandle(gmmFig)
        sumGamma = sum(gammas) / numVoxels;
        for i = 1:numGaussians
            gaussianEvals(:, i) =  (1 / sqrt(2*pi*variances(i))) * exp(-0.5 * (binCenters - means(i)).^2 / variances(i)) * sumGamma(i) * binWidth;
        end
        figure(gmmFig)
        clf
        grid
        bar(binCenters, histogram);
        hold on
        plot(binCenters, gaussianEvals', 'g');
        plot(binCenters, sum(gaussianEvals, 2), 'r', 'lineWidth', 1);
        hold off
        drawnow
    end
    
    for i = 1:numGaussians
        squaredDists(:, i) = (data - means(i)).^2;
        likelihoods(:, i) = (1 / sqrt(2*pi*variances(i))) * exp(-0.5 * squaredDists(:,i) / variances(i)) .* gammas(:, i);
    end
    
    normalizer = sum(likelihoods, 2);
    
    minLogL = -sum(log(normalizer + eps));
    costHist(end + 1) = minLogL;%#ok
    
    if(costHist(end) > costHist(end-1))
        warning('flaf');
    end
    
    relChange = (costHist(end - 1) - costHist(end)) / numVoxels;
    if abs(relChange) < convCriterion
        if ~isempty(strfind(verbosity, 'debug'))
            fprintf('\tGGM fit converged in %i iterations\n', mmIter);
        end
        converged = 1;
        break
    end
    
    if ishandle(costFig)
        figure(costFig)
        plot(costHist(2:end));
        drawnow
    end
    
    posteriors = likelihoods ./ repmat(normalizer + eps, [1, numGaussians]); 
    posteriorsSum = sum(posteriors);
    
    % classes we trust, to avoid collapses
    trusted = find(posteriorsSum >  (1e-3 * sum(posteriorsSum)));
    
    if(~all(trusted))
        warning('Classes are collapsing');
    end
    % update weights
    weights = posteriorsSum / sum(posteriorsSum);
    
    if(abs(sum(weights) - 1) > 1e-7)
        error('Weights do not sum up to one')
    end
           
    if updateMeans 
        if(equidistantMeans)
            lhs = zeros(2, 2);
            lhs(1, 1) = sum(sum(posteriors) .* (alphas').^2);
            lhs(2, 1) = sum(sum(posteriors) .* (alphas') .* (1-alphas'));
            lhs(1, 2) = lhs(2, 1);
            lhs(2, 2) = sum(sum(posteriors) .* (1-alphas').^2);
            
            rhs = [(data' * posteriors) * alphas; (data' * posteriors) * (1-alphas)];
            
            independentMeans = lhs \ rhs;
            means = ([alphas, (1-alphas)] *  independentMeans)';
        else
            newmeans = sum(posteriors .* repmat(data,[1, numGaussians])) ./ (posteriorsSum + eps);
            means(trusted) = newmeans(trusted);
        end
    end
    
    if updateVariance
        for i = 1:numGaussians
            squaredDists(:, i) = (data - means(i)).^2;
        end
        if(sameVariance)
            variance = posteriors(:)' * squaredDists(:) / sum(posteriorsSum);
            variances = repmat(variance, [1, numGaussians]);
        else
            newvariances = sum(posteriors .* squaredDists) ./ (posteriorsSum + eps);
            newvariances = max(newvariances, 0.001);
            variances(trusted)= newvariances(trusted);
        end
    end
    
    % update gammas and recompute posteriors taking into account the MRF property. We only do one MRF update to avoid local minima
    if betaMRF
        [gammas, ~] = updateWithMRF(data, posteriors, means, variances, weights, betaMRF, updatePatternMRF, maxIters, convCriterion, verbosity);
    end
    
end

if (~converged && ~isempty(strfind(verbosity, 'debug')))
    disp(['GMM stopped: Maximum number (', num2str(maxIters), ') of iterations reached'])
end

if(nargout > 4)
    Sik = posteriors ./ repmat(variances, [numVoxels, 1]);
    Si = sum(Sik, 2);
end

if(nargout > 5)
    expectation = sum(Sik .* repmat(means, [numVoxels, 1]), 2) ./ (Si + eps);
end

end

function [gammas, posteriors] = updateWithMRF(data, posteriors, means, variances, weights, beta, updatePatternMRF, maxIters, convCriterion, verbosity)

numClasses = size(posteriors, 2);
mask = updatePatternMRF{1};
numVoxels = sum(mask(:) > 0);

nonPadded = zeros(size(mask));
padded = zeros(size(nonPadded) + 2);
gammas = zeros(size(posteriors));

likelihoods = zeros(numVoxels, numClasses);

costHist = 1 / eps;

for iterMRF = 1:1
    for i = 2:3
        idx = updatePatternMRF{i};
        for j = 1:numClasses
            % Compute gamma eq (4.23)
            nonPadded(mask) = 1 - posteriors(:, j);
            padded(2:end-1, 2:end-1, 2:end-1) = nonPadded;
            
            posteriorSumNeighborhood = ...
                padded((2:end-1) - 1, (2:end-1), (2:end-1)) + ...
                padded((2:end-1), (2:end-1) - 1, (2:end-1)) + ...
                padded((2:end-1), (2:end-1), (2:end-1) - 1) + ...
                padded((2:end-1) + 1, (2:end-1), (2:end-1)) + ...
                padded((2:end-1), (2:end-1) + 1, (2:end-1)) + ...
                padded((2:end-1), (2:end-1), (2:end-1) + 1);
            
            % MRF priors (gamma), Equation (4.23)
            gammas(:, j) = log(weights(j)) - beta * posteriorSumNeighborhood(mask);
            
        end
        
        % compute gammas
        gammaExp = exp(gammas);
        gammas = gammaExp ./ repmat(sum(gammaExp, 2) + eps, [1, numClasses]);
       
        %Re-estimate posteriors and normalize
        for j = 1:numClasses
            likelihoods(idx, j) = 1/sqrt(2 * pi * variances(j)) * exp(-0.5 * (data(idx) - means(j)).^2 / variances(j)).* gammas(idx, j);
        end           
        
        normalizer = sum(likelihoods, 2);
        posteriors = likelihoods ./ repmat(normalizer + eps, [1, numClasses]); 
    end
    
%     costHist(end+1) = -sum(log(normalizer + eps)); %#ok
%     relChange = (costHist(end - 1) - costHist(end)) / numVoxels;
%     if abs(relChange) < convCriterion
%         if ~isempty(strfind(verbosity, 'debug'))
%             fprintf('\tMRF gamma estimation converged in %i iterations\n', iterMRF);
%         end
%         break
%     end
end
        


end
