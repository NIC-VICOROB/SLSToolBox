function [means, variances, weights, minLogL, posteriors, Si, expectation] = fitGMMwithPrior(...,
    data, ...
    prior, ...
    priorMap, ...
    means, ...
    variances,...
    weights, ...
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
if nargin < 7
    maxIters = 1000;
end
if nargin < 8
    convCriterion = 1e-5;
end
if nargin < 9
    trustWeights = 1e-3;
end
if nargin < 10
    trustVariance = 1e-3;
end
if nargin < 11
    numBins = 1024;
end
if nargin < 12
    verbosity = 'none';
end
if nargin < 13
    gmmFig = -1;
end
if nargin < 14
    costFig = -1;
end

% initialize
numGaussians = length(means);
numClasses = size(prior, 2);
numVoxels = numel(data);
posteriors = zeros(numVoxels, numGaussians);

% history of minLogLikelihoods - used to determine convergenge between iterations
costHist = 1/eps;
converged = false;

if ishandle(gmmFig)
    gaussianEvals = zeros(numBins, numGaussians);
    classEvals = zeros(numBins, numClasses);
    
    [histogram, binCenters] = hist(data, numBins);
    histogram = histogram / sum(histogram);
    binWidth = binCenters(2) - binCenters(1);
    
    priorFrequency = sum(prior) / numVoxels;
end

% Performing EM
for mmIter = 1:maxIters
    
    % plot the state of affairs
    if ishandle(gmmFig)
        figure(gmmFig)
        for k = 1:numGaussians
            gaussianEvals(:, k) =  (1 / sqrt(2*pi*variances(k))) * ...
                exp(-0.5 * (binCenters - means(k)).^2 / variances(k)) * weights(k) * priorFrequency(priorMap(k));
        end
        gaussianEvals = gaussianEvals * binWidth;
        for k = 1:numClasses
            classEvals(:,k) = sum(gaussianEvals(:,k == priorMap),2);
        end
        clf
        grid
        bar(binCenters, histogram);
        hold on
        plot(binCenters, classEvals', 'lineWidth', 1);
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
        posteriors(:, k) =  (1 / sqrt(2*pi*variances(k))) * exp(-0.5 * ((data - means(k)).^2) / variances(k)) ...
                            .* prior(:,priorMap(k)) * weights(k);
    end
    clear imageModel
  %  likelihoods = likelihoods .* prior(:,priorMap);
    
    normalizer = sum(posteriors, 2);

    for i = 1:numGaussians
        posteriors(:,i) = posteriors(:,i) ./ (normalizer + eps);
    end
    
    minLogL = -sum(log(normalizer + eps));
    costHist(end + 1) = minLogL;%#ok
    
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
        axis tight
        title('Cost');
        xlabel('Iteration');
        ylabel('Cost');
        drawnow
    end
   
    posteriorsSum = sum(posteriors);
    posteriorsSumSum = zeros(numClasses, 1);
    
    for g = 1:numGaussians
        l = priorMap(g);
        posteriorsSumSum(l) = posteriorsSumSum(l) + sum(posteriorsSum(g));
    end
    
    % update weights
    for g = 1:numGaussians
        l = priorMap(g);
        weights(g) = posteriorsSum(g) / posteriorsSumSum(l);
    end
    
    % verify that \sum w_gl = 1
    for l = 1:numClasses
        if(abs(sum(weights(priorMap == l)) - 1) > 1e-7)
            error('Weights do not sum up to one')
        end
    end
    
    % don't update gaussians that contribute less than trustCriterion to the posterior ( which should == 1)
    trusted = find(weights >  trustWeights);
    if(~all(trusted))
        warning('Classes are collapsing');
    end
    % update means
    % newmeans = sum(posteriors .* repmat(data,[1, numGaussians])) ./ (posteriorsSum + eps);
    newmeans = means;
    %memory consideration
    for k = 1:numGaussians
        newmeans(k) = (posteriors(:,k)' * data) ./ (posteriorsSum(k) + eps);
    end
    means(trusted) = newmeans(trusted);
    
    % update variances
%     newvariances = sum(posteriors .* squaredDists) ./ (posteriorsSum + eps);
%     newvariances = max(newvariances, trustVariance);

    % memorys
    newvariances = variances;   
    for k = 1:numGaussians
        squaredDists = ((data - means(k)).^2);
        newvariances(k) = (posteriors(:,k)' * squaredDists) ./ (posteriorsSum(k) + eps);
        newvariances(k) = max(newvariances(k), trustVariance);  
    end   
    variances(trusted)= newvariances(trusted);
end

if (~converged && ~isempty(strfind(verbosity, 'debug')))
    disp(['GMM stopped: Maximum number (', num2str(maxIters), ') of iterations reached'])
end

if(nargout > 5)
    Si = zeros(numVoxels, 1);
    expectation = zeros(numVoxels, 1);
    for k = 1:numGaussians
        Sik = posteriors(:,k) / variances(k);
        
        Si = Si + Sik;
        expectation = expectation + Sik * means(k);
    end
    expectation = expectation ./ Si;
end

end