function threshold = biModalThreshold(data, useLogTransform, method, numberOfBins, verbosity, thresFig)
%

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


if nargin < 2
    useLogTransform = false;
end
if nargin < 3
    method = 'mixture';
end
if nargin < 4
    numberOfBins = 1024;
end
if nargin < 5
    verbosity = 0;
end
if nargin < 6
    thresFig = -1;
end


% Get histogram in log domain
mask = (data > 0);
if useLogTransform
    data = log(data(mask));
else
    data = data(mask);
end

range = [min(data(:)), max(data(:))];
rangeSpan = range(2) - range(1);
binWidth = rangeSpan / (numberOfBins-1);
binCenters = range(1) + (0:numberOfBins-1) * binWidth;

histogram = hist(data, binCenters);

method = lower(method);
        

switch method
    case 'mixture'
        % Fit a two-component mixture model to the data
        
        equilateralMeans = false;
        updateMeans = true;
        sameVariance = false;
        updateVariance = true;
        maxIters = 500;
        convCriterion = 1e-6;
        numBins = numberOfBins;
        
        % Remove intensity outliers when initialising
        percentile = 0.05;
        tmp = cumsum(histogram / sum(histogram));
        firstRealBin = find(tmp > percentile, 1);
        lastRealBin = find(tmp < (1 - percentile) , 1, 'last');
        
        range = [binCenters(firstRealBin), binCenters(lastRealBin)];
        rangeSpan = range(2) - range(1);
        
        % 2 Gaussians, since we do background foreground
        numGaussians = 2;

        gaussianWidth = rangeSpan / numGaussians;
        % Initialize means based on min/max of data
        means = range(1) + 0.5 * gaussianWidth + (0:numGaussians-1) * gaussianWidth;
        variance = gaussianWidth^2;
        variances = repmat(variance, [1, numGaussians]);
        weights = ones(1, numGaussians) / numGaussians;
        
        numVoxels = length(data);
        
       % intensities = data(data >= firstRealBin & data <= lastRealBin);
        [means, variances, weights] = fitGMM(...
            data, ...
            means, ...
            variances, ...
            weights, ...
            equilateralMeans, ...
            updateMeans, ...
            sameVariance, ...
            updateVariance, ...
            maxIters, ...
            convCriterion, ...
            numBins, ...
            verbosity, ...
            thresFig);
        
        % Now get the binNumber that seperates the two
        % Pre-calculate mahalanobis distances (useful for EM update equation of the variance)
        squaredDistances = zeros(numVoxels, numGaussians);
        likelihoods = zeros(numVoxels, numGaussians);
        
        for gaussianNumber = 1:numGaussians
            variance = variances(gaussianNumber);
            squaredDistances(:, gaussianNumber) = (data - means(gaussianNumber)).^2;
            likelihoods(:, gaussianNumber) = exp(-0.5 * squaredDistances( : , gaussianNumber) / variance) / sqrt(2 * pi * variance);
        end
        
        if ishandle(thresFig)
            histogram = histogram / sum(histogram);
            binWidth = binCenters(2) - binCenters(1);
            gaussianEvals = zeros(numBins, numGaussians);
            for i = 1:numGaussians
                gaussianEvals(:, i) =  (1 / sqrt(2*pi*variances(i))) * exp(-0.5 * (binCenters - means(i)).^2 / variances(i)) * weights(i) * binWidth;
            end
            figure(thresFig)
            clf
            grid
            bar(binCenters, histogram);
            hold on
            plot(binCenters, gaussianEvals', 'g');
            plot(binCenters, sum(gaussianEvals, 2), 'r', 'lineWidth', 2);
            hold off
        end
        
        % Calculate posteriors
        posteriors = likelihoods .* repmat(weights, [numVoxels, 1]);
        normalizer = sum(posteriors, 2) + eps;
        posteriors = posteriors ./ repmat(normalizer, [1, numGaussians]);
        background = posteriors(:, 1) - posteriors(:, 2);
        
        % figure
        %  plot( binCenters, posteriors( :, 1 ) )
        %  hold on
        %  plot( binCenters, posteriors( :, 2 ), 'r' )
        %  figure
        %  plot( binCenters, background );
        
        % we need to identify the best intersection between the two Gaussians; if we do not do this, we risk that the
        % method incorrectly thresholds at max intensities (very far from the mean of both distributions
        background = background(data >= means(1) & data <= means(2));
        data = data(data >= means(1) & data <= means(2));
        
        % Determine threshold
        [~, minIdx] = min(abs(background));   
        threshold = data(minIdx);
        if useLogTransform
            threshold = exp(data(minIdx));
        end
        
    case 'otzu'
        % Exhaustive search with Otzu method
        bestThreshold = -Inf;
        bestInterClassVariance = -Inf;
        for threshold = binCenters
            histogramOfLeft = histogram .* (binCenters < threshold);
            histogramOfRight = histogram .* (binCenters >= threshold);
            meanOfLeft = binCenters * histogramOfLeft' / sum(histogramOfLeft + eps);
            meanOfRight = binCenters * histogramOfRight' / sum(histogramOfRight + eps);
            weightOfLeft = sum(histogramOfLeft) / sum(histogram);
            weightOfRight = sum(histogramOfRight) / sum(histogram);
            
            interClassVariance = weightOfLeft * weightOfRight * (meanOfLeft - meanOfRight)^2;
            if interClassVariance > bestInterClassVariance
                bestInterClassVariance = interClassVariance;
                bestThreshold = threshold;
            end
            
        end 
        
        threshold = bestThreshold;    
        
        % plot the thresholded figure
        if ishandle(thresFig)
            figure(thresFig)
            histogramOfLeft = histogram .* (binCenters < threshold);
            histogramOfRight = histogram .* (binCenters >= threshold);
            bar(binCenters,  histogramOfLeft, 'facecolor', 'b', 'edgecolor', 'b')
            hold on
            bar(binCenters,  histogramOfRight, 'facecolor', 'r', 'edgecolor', 'r')
        end
        
        if useLogTransform
            threshold = exp(bestThreshold);
        end     
    case 'n3'
        % I guess this is quite similar to Otzu - but in order to ensure that we do things the same way, here follows
        % yet another implementation that does things the N3 way.
        
        N = sum(histogram);
        %meanValue = sum(binCenters .* histogram) / N;
        meanValue = sum((binCenters + 0.5*binWidth) .* histogram) / N;
        varMax = 0;
        i = 1;       
        offset = binWidth * 0.5;
        %offset = 0;

        zeroMoment0  = histogram(1)/N;
        firstMoment0 = (binCenters(1) + offset) * histogram(1)/N;

        for k = 2:numberOfBins
            zeroMoment1  = zeroMoment0 + histogram(k)/N;
            firstMoment1 = firstMoment0 + (binCenters(k) + offset)*histogram(k)/N;
            
            if zeroMoment1 > 0 && zeroMoment1 < 1
                var = (meanValue*zeroMoment1 - firstMoment1)^2 / (zeroMoment1*(1 - zeroMoment1));
                if var > varMax
                    i = k;
                    varMax = var;
                end
            end
            
            zeroMoment0 = zeroMoment1;
            firstMoment0 = firstMoment1;
        end
        
        threshold = binCenters(i) + offset;  
        
         % plot the thresholded figure
        if ishandle(thresFig)
            figure(thresFig)
            histogramOfLeft = histogram .* (binCenters < threshold);
            histogramOfRight = histogram .* (binCenters >= threshold);
            bar(binCenters,  histogramOfLeft, 'facecolor', 'b', 'edgecolor', 'b')
            hold on
            bar(binCenters,  histogramOfRight, 'facecolor', 'r', 'edgecolor', 'r')
        end       
        
        if useLogTransform
            threshold = exp(threshold);
        end
        
        % Conform with N3 argument parsing
       % threshold = floor(threshold);    
    otherwise
        error('threshold method undefined');
end

end








