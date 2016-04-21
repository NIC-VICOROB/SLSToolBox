function [means, variances, weights, minLogL, dataEstimate, elapsedTimeWeights] = fitGMMWiener(...,
    data, ...
    numGaussians, ...
    fwhm, ...
    wienerFiltNoise, ...
    numBins, ...
    useInterpHist, ...
    useInterpTable, ...
    doPad, ...
    gmmFig)
% Description: Fits a Gaussian mixture model to data by means of N3 histogram sharpening (Wiener Filtering)
%
% Call: [means, variances, weights, minLogL, voxelWeights, dataEstimate] = fitGMM(data, numGaussians, fwhm,
%           wienerFiltNoise, numBins, figToPlotTo)
%
% Input:    data:               Data to fit the GMM to
%           numGaussians:       Number of Gaussians that should be used for the histogram sharpening (number of bins)
%           fwhm:               'Full Width at Half Maximum' - parameter that is used to compute Gaussian variance
%           wienerFiltNoise:    Noise value used in the Wiener Filtering
%           numBins:            bins to plot histogram with
%           useInterpHist:      Create weighed histogram by means of interpolation
%           useInterpTable:     Lookup voxel expectations by means of histogram interpolation
%           doPad:              Pad the gaussians to conform with N3 wiener deconvolution in Fourier domain
%           figToPlotTo:        Show the Mixture Model fit after the histogram has been sharpened Wiener Filtering style
%
% Output:   means:              Gaussian means that optimize likelihood of data
%           variances:          Gaussian variances that optimize likelihood of data
%           weights:            Gaussian weights that optimize likelihood of data
%           minLogL:            The negative loglikelihood at convergence
%           voxelweights:       The weights of all voxels in fit (posterior probability over variance)
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

if nargin < 3
    fwhm = 0.15;
end
if nargin < 4
    wienerFiltNoise = 0.01;
end
if nargin < 5
    numBins = 1024;
end
if nargin < 6
    useInterpHist = 1;
end
if nargin < 7
    useInterpTable = 1;
end
if nargin < 8
    doPad = 1;
end
if nargin < 9
    gmmFig = -1;
end

%doPad = 0;
numVoxels = length(data);

% Define range of intensities that our non-zero data spans from DHistogram.cc
range = [min(data(:)), max(data(:))];
rangeSpan = range(2) - range(1);
binWidth = rangeSpan / (numGaussians - 1);

% AHA! N3 DOES use the parzen window histogram forming, which depends on WHistogram, NOT DHistogram
range(1) = range(1) - 0.5*binWidth;
range(2) = range(2) + 0.5*binWidth;
binEdges = range(1):binWidth:range(2);
binCenters = range(1) + 0.5 * binWidth + (0:(numGaussians-1)) * binWidth;

% Form histogram from WHistogram.h
if useInterpHist
    histogram = interpHist(data, numGaussians);
    if ~useInterpTable
        % For looking up - same as combining what they do in sharpen_hist.cc + minclookup
        [~, binNumbers] = histc(data, binEdges);
    end
else
    % Form histogram given previously computed bin edges
    [histogram, binNumbers] = histc(data, binEdges);
    
    % histc creates one bin too many with the extreme maxima - therefore we need to cut the histogram by one
    histogram(end-1) = histogram(end-1) + histogram(end);
    histogram = histogram(1:end-1);
end

if ~useInterpTable
    % and update the intensities with the new bin
    binNumbers(binNumbers == numGaussians+1) = numGaussians;

    % % This should NEVER happen with the way we handle the histogram generation
    if(any(binNumbers == 0))
        error('An intensity was not properly included in histogram');
    end
end

% slope as defined in sharpen_hist.cc - affects variances
% binCenters and means are also updated since histogram has changed and histogram needs padding
slope = (max(binCenters(:)) - min(binCenters(:))) / (numGaussians - 1);
variance = 1 / (8 * log(2)) * (fwhm / slope)^2* binWidth^2;
variances = repmat(variance, [1, numGaussians]);

% to mimic the fourier operation, we pad and expand
%binCenters = range(1) + 0.5 * binWidth + ((-offset):(numBins-1+offset)) * binWidth;
means = binCenters;
numGaussiansPadded = numGaussians;
if doPad
    numGaussiansPadded = pow2(ceil(log(numGaussians)/log(2.0))+1);
    offset = 0.5 * (numGaussiansPadded - numGaussians);
    oldHistogram = histogram;
    histogram = zeros(numGaussiansPadded, 1);
    histogram((offset+1):(numGaussians+offset)) = oldHistogram;
    binCenters = range(1) + 0.5 * binWidth + ((-offset):(numGaussians-1+offset)) * binWidth;
    means = binCenters;
    variances = repmat(variance, [1, numGaussiansPadded]);
end

% Calculate basis functions
squaredDists = zeros(numGaussiansPadded, numGaussiansPadded);
likelihoods = zeros(numGaussiansPadded, numGaussiansPadded);

% fwhm = fwhm / slope;
% factor = 4*log(2)/(fwhm*fwhm);
% scale = 2*sqrt(log(2)/pi) / fwhm;

for i = 1:numGaussiansPadded
    %     ds = ((1:200) - i).^2;
    %likelihoods(:, i) = scale*exp(-ds*factor);
    squaredDists(:, i) = (binCenters'-means(i)).^2;
    likelihoods(:, i) = exp(-0.5 * squaredDists(:, i) / variances(i)) / sqrt(2 * pi * variances(i));
end

% Perform least-squares fitting

% This step is necessary to mimic N3 because we adjusted variance according the binWidth
% (which we did to associate bin centers with our log-tranformed intenties, whereas N3 uses a simple intensity decoupled loop counter)
likelihoods = likelihoods * binWidth;

lhs = likelihoods' * likelihoods + diag(wienerFiltNoise*ones(numGaussiansPadded, 1));
rhs = likelihoods' * histogram;

startTimeWeights = tic;
weights = lhs \ rhs;
elapsedTimeWeights = toc(startTimeWeights);
% LS fitting yields negative weights, which is obviously wrong in a Gaussian mixture
% N3 hacks this by zeroing negative weights
weights(weights < 0) = 0;

weights = weights' / sum(weights);

% Calculate minLogLikelihood
posteriors = likelihoods .* repmat(weights, [numGaussiansPadded, 1]);

normalizer = sum(posteriors, 2);
posteriors = posteriors ./ repmat(normalizer + eps, [1, numGaussiansPadded]);

%Construct a lookup table of where every intensity wants to go as per EM
if nargout > 4
    Sik = posteriors ./ repmat(variances, [numGaussiansPadded, 1]);
    Si = sum(Sik, 2);
    
    predictedImageLookup = sum(Sik .* repmat(means, [numGaussiansPadded, 1]), 2) ./ (Si + eps);
    if doPad
        predictedImageLookup = predictedImageLookup((offset+1):(numGaussians+offset));
    end
    if useInterpTable
        dataEstimate = interpTable(data, predictedImageLookup);
    else
        dataEstimate = predictedImageLookup(binNumbers);
    end
end

if doPad
	means = means((offset+1):(numGaussians+offset));  
    variances = variances((offset+1):(numGaussians+offset));    
    weights = weights((offset+1):(numGaussians+offset));    
end


% compute cost from voxels
squaredDists = zeros(numVoxels, numGaussians);
likelihoods = zeros(numVoxels, numGaussians);
for i = 1:numGaussians
    squaredDists(:, i) = (data - means(i)).^2;
    likelihoods(:, i) = (1 / sqrt(2*pi*variances(i))) * exp(-0.5 * squaredDists(:,i) / variances(i)) * weights(i);
end
normalizer = sum(likelihoods, 2);

minLogL = -sum(log(normalizer + eps));


% plot stuff
if ishandle(gmmFig)
    [histogram, binCenters] = hist(data, numBins);
    histogram = histogram / sum(histogram(:));
    binWidth = binCenters(2) - binCenters(1);
    
    gaussPoints = linspace(min(data(:)), max(data(:)),512);
    gaussianEvals = zeros(512, numGaussians);
    for i = 1:numGaussians
        gaussianEvals(:,i) =  (1 / sqrt(2*pi*variances(i))) * exp(-0.5 * (means(i) - gaussPoints).^2 / variances(i)) * weights(i) * binWidth;
    end
    figure(gmmFig)
    clf
    grid
    bar(binCenters, histogram);
    hold on
    plot(gaussPoints, gaussianEvals', 'g','lineWidth', 1);
    plot(gaussPoints, sum(gaussianEvals, 2), 'r', 'lineWidth', 2);
    xlabel('Bins')
    xlim([1,7]);
    title('Log-intensity histogram');
    ylabel('Normalized Voxel Count');
    set(gca, 'fontsize', 6)
    axis tight
    hold off
    drawnow
end

end

function histogram = interpHist(data, numBins)

numVoxels = numel(data);

% Define range of intensities that our non-zero data spans from DHistogram.cc
range = [min(data(:)), max(data(:))];
rangeSpan = range(2) - range(1);
binWidth = rangeSpan / (numBins - 1);
range(1) = range(1) - 0.5*binWidth;
range(2) = range(2) + 0.5*binWidth;
rangeSpan = range(2) - range(1);

% AHA! N3 DOES use the parzen window histogram forming, which depends on WHistogram, NOT DHistogram
% histogram = zeros(numBins, 1);
histogram = zeros(numBins, 1);

%
loc = numBins * (data - range(1)) / rangeSpan + 1;
index = floor(loc);
offset = loc - index - 0.5;

for i = 1:numVoxels
    if (abs(offset(i)) < 1e-7)
        histogram(index(i)) = histogram(index(i)) + 1;
    elseif (offset(i) > 0 && index(i) < numBins)
        histogram(index(i)) = histogram(index(i)) + 1 - offset(i);
        histogram(index(i) + 1) = histogram(index(i) + 1) + offset(i);
    elseif (index(i) > 1)
        histogram(index(i)) = histogram(index(i)) + 1 + offset(i);
        histogram(index(i)-1) = histogram(index(i)-1) - offset(i);
    else     
        error('this should not happen');
    end
end

if abs(sum(histogram) - numVoxels) > 1e-7
    error('histogram does not count up to the number of voxels accurately enough');
end

end


function data = interpTable(data, tableValues)

numVoxels = numel(data);

% % sort data to guarantee the following works
[data, idx] = sort(data);

% 
% % 0 -> 1 mappings
dataRange = [min(data), max(data)];
dataRangeSpan = dataRange(2) - dataRange(1);
index = (data - dataRange(1)) / dataRangeSpan;

% this is from sharpen_hist.cc
table = (0:(numel(tableValues)-1)) / (numel(tableValues)-1);

% handle cases where the value is too far to the left relative to the table
mask = index < table(1);
data(mask) = tableValues(1);

i = sum(mask) + 1;
j = 2;

% now look up values that are within the table limits
while table(j) < table(end)
    % normalizer, given table values
    oneOverSpan = 1 / (table(j) - table(j-1));
    
    % look up values until we need to move forward in the table again
    while index(i) < table(j)
        % interpolation fraction between two table points
        fraction = (index(i) - table(j-1)) * oneOverSpan;
        
        % interpolated data
        data(i) = (1-fraction)*tableValues(j-1) + fraction*tableValues(j);
        
        % next data value
        i = i + 1;
    end
    
    % move forward in table
    j = j + 1;
end

% handle cases where the value is too far to the right relative to table
data(i:numVoxels) = tableValues(end);

% reorder data to original order
data(idx) = data;

end

