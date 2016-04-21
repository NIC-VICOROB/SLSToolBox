function [correctedData, estimatedBiasField, mask, infoStruct, segmented] = iic(data, varargin)
% Description: Estimates a multiplicative intensitity inhomogeneity field from data by means of Bayesian inference
%
% CALL: [correctedData, biasField, mask, info] = iic(data, OPTIONAL)
%
% INPUT:
%   data:                       MR image (3D matrix) to be corrected. If the matrix is 4D, the fourth dimension contains timepoints.
%
% OPTIONAL INPUT:
%   'stepSize':                 Step size of data (voxel width, height and depth) in mm. Default: [1,1,1].
%
%   'desiredStepSize':          The step size you want to estimate at (used for downsampling the input data) in mm. Default: 4.
%
%   'backgroundThresholdType':  Method used to compute background treshold. Used to mask out low intensity voxels.
%                               Available options: 'none', 'mixture', 'otzu' (similar to N3), 'superVoxel'.
%                               If no mask and no prior is provided, the method defaults to 'otzu', else 'none'.
%                               'mixture' is a standard 2 component gaussian mixture model fit.
%                               'otzu' is a speedy approximation to 'mixture', which works really well in practice.
%                               'superVoxel': See mixture fit details. After fitting supervoxels to the image, 'otzu' is used to segment 
%                               supervoxel means into foreground and background,Background is then filtered away, and foreground is used for
%                               supervoxel mixture fit initialization.
%
%   'maxNumberIterations':      Maximum number of iterations that the method is allowed to run. Default: 1000. Original N3 default: 50.
%
%   'maxNumberGMMIterations':   Maximum number of iterations that the mixture fitting component is allowed to run. Default: 100.
%
%   'mask':                     A mask specifying voxels to be used for intensity inhomogeneity estimation.
%                               Must be same dimensions as data.
%
%   'prior':                    An IxP sized matrix of P spatially smooth priors on tissue classes or structures for I image voxels.
%                               All voxel priors are checked (and normalized if \sum_p v_p^i > 1) such that \sum_p v_p^i <= 1.
%                               A garbage/background class is computed from the priors to ensure that \sum_p v_p^i = 1.
%
%   'priorThreshhold':          Percentage threshhold for masking out background voxels when a prior is provided. Default: 0.99.
%                               This mask is combined with a thresholding of zero-intensity voxels.
%
%   'superVoxelStepSize'        The desired spacing in mm between supervoxels. The number of supervoxels is derived from this + the data.
%                               The result is P (prod(superVoxels)) Gaussians that are fitted, one each, to P super voxels with centers
%                               that are initially equidistantly distributed along each direction.
%                               Super voxel centers that falls outside the mask are filtered away.
%                               if superVoxelStepSize is specified, numGaussians becomes  = P, regardless of input (for now...)
%                               Default: 0 (disabled).
%
%   'superVoxelTheta'           Weight (in %) between intensity and spatial variance. Only applies for sameVariance. Default is 100%.
%
%
%   'numberOfGaussians':        Number of Gaussians to use in mixture component.
%                               If a prior is specified, numberOfGaussians is a vector containing the number of Gaussians to use per class.
%                               Default: 3 gaussians in a single class, original N3 uses a default value of 200 (hidden).
%
%   'equidistantMeans':         Default: false. If true, Gaussian means are restricted to equidistant placement.
%                               Original N3: true.
%
%   'updateMeans':              Default: true. If set to false, means are never updated after initialization.
%
%   'sameVariance':             Default: false. If set to true, Gaussian variances are restricted to equal value.
%                               Original N3 default: true.
%
%   'updateVariance':           Default: true. If set to false, variances are never updated after initialization
%                               
%                               Original N3: false.
%
%   'varianceInit':             Default: 0 (disabled). Specify starting value for Gaussian variance. 
%                               Defaults to squared distance between means if not specified.
%
%   'gmmFitType':               Method used to fit the mixture model.
%                               Options: 'mixture' (default), 'mixtureWithLabelPrior', mixtureWithSuperVoxel', 'wiener'.
%                               Original N3 uses 'wiener', a heuristic gaussian mixture fitting method performed by a
%                               regularized least-squares method (Wiener deconvolution).
%
%   'convergenceType':          Correction method convergence type. Options: 'cost' (default), uses change in
%                               cost per voxel relative to cost starting value (computed from initialized parameters).
%                               N3 fitting implies 'stdDifferenceField', which uses standard deviation
%                               of difference in field estimates between current and previous iteration.
%
%   'convergence':              Default: 1e-5. Convergence threshold. Original N3 uses a default value of 1e-3.
%
%   'convergenceGMM':           Default: 1e-5. Convergence threshold for Gaussian mixture. Only applies if 'gmmFitType': 'mixture'.
%
%   'smoothingType':            Type of smoothing used on intensity residual, in order to conform with the assumption of a smooth,
%                               multiplicative bias. Options: 'cosine', 'N3' (default). 'N3' is really a bspline scheme with
%                               equidistant control points and regularization on the partial second order derivatives (bending
%                               energy).
%
%   'smoothingDistance':        Default: 50. Distance in mm which controls how flexible the bias is. Implicitly controls the number
%                               of basis functions used for smoothing.
%                               Note that this parameter can be specified per dimension (e.g., [50, 60, 70]), which allows (discrete)
%                               control with how flexible the bias is along each dimension.
%
%   'smoothingRegularization':  Default value: 1e-5. A scalar that controls how heavily we regularize (penalize curvature).
%                               Must be specified when using bspline smoothing ('N3' default).
%                               Original N3 uses 1e-7.
%
%   'voxelScaleRegularization': Default: true. Scale the regularization by the number of voxels (balances reg term vs dataterm).
%                               Original N3 fitting uses a default value of 'true'.
%
%   'smoothExponentiatedBias':  (Disabled). A smoothing performed by the original N3 method of the final
%                               exponentiated bias field estimate. Generally found to do very little, and it might hurt
%                               performance (e.g., CJV between WM and GM). Default in 'strict' N3, otherwise disabled.
%
%   'fwhm':                     Default: 0.15. N3 parameter used to control variance of Gaussians in mixture. Implies
%                               'gmmFitType': 'wiener'.
%
%   'wienerFilterNoise':        Default: 0.1. N3 hidden parameter used to control noise used by wiener fitting of
%                               gaussian mixture.
%
%   'mimicN3':                  Options: 'no' (default), 'model', 'strict'.
%                               'model' implies N3 mimicking and enables all model essential N3 parameters.
%                               'strict' further mimics N3 by enabling various tweaks and hacks.
%                               Parameters that can be user controlled in the original N3 are NOT set to the original defaults.
%                               Note that you can configure N3 like behavior by setting other flags yourself.
%
%   'histogramBins':            Default: 256. Number of bins used when plotting histogram and mixture fitting.
%
%   'showPlots':                Default: false. True enables a number of plots that are useful to debug and inspect
%                               progress of correction.
%
%   'verbosity':                Options: 'none' (default), 'info' and 'debug'. Print relevant information to commandline
%                               during progress of correction.
%
% OUTPUT:
%   correctedData:              Intensity homogeneity corrected data at original resolution.
%
%   biasField:                  The estimated intensitiy inhomogeneity field at original resolution.
%
%   mask:                       The mask of voxels that are used for bias estimation at original resolution. Useful if
%                               mask is computed using background thresholding.
%
%   info:                       A struct containing various history vectors, useful for debugging and inspection of
%                               method behavior.
%
% Author:       Christian Thode Larsen (cthla@dtu.dk), christian.thode.larsen@gmail.com
% Literature:   Leemput et al 1999. Automated Model-based Bias Field Correction of MR Images of the Brain.
%               Larsen et al 2014. N3 Bias Field Correction Explained as a Bayesian Modeling Method.
%               Larsen 2015. Development and Application of Tools for MRI Analysis. PhD Thesis.

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

% EDIT 8/2/2013
% Changed cost convergence criterion to allow wiener filter convergence with abs(relativeCost)
% Modified to operate on proper number of Gaussians for mimicN3: 200
%
% EDIT 9/17/2013
% Introduced N3 b-spline smoothing scheme to have a better 1:1 mapping between mimicN3 and the original method
% Minor edits to how the masks are generated; data is clamped to [1; 1.7e308], and a background mask is identified for
% values of 1
%
% EDIT 9/19/2013
% Changed how histogram edges and centers are computed to mimic N3. Seems more precise and compact, and updated in every
% iteration (before safety margins were used and the histogram bin edges and centers were never updated)
% MILESTONE: Performance for mimicN3 is now very similar to original N3
%
% EDIT 9/21/2013
% fixed N3 slope computation and moved within loop, since it depends on histogram max and min bin centers - which change in
% every iteration. This affects the placement of means.
% Introduced padding of the histogram according to sharpen_hist.cc (200 bins are zero-padded to 512).
%
% EDIT 11/14/2013
% Total code clean-up, several bugs eliminated
% Cost is now computed in a more robust and accurate manner + cost of regularizer supported
% Supports N3 scaling of lambda (their code multiplies by numVoxels and not 2*variance)
%
% EDIT 11/21/2013
% Implemented N3 convergence scheme (std dev, and NOT coefficient of variation as claimed in article. This yields approx
% same number of iters as N3
%
% EDIT 04/01/2014
% Full implementation of bsplines + regularization the N3 way. N3 I/O dependency removed - Milestone for MICCAI 2014
%
% EDIT 07/17/2014
% Too many edits to summarize everything... FAST smoothing using anisotropic basis functions supported (yay!).
%
% EDIT 07/21/2014
% Input/output cleanup
%
% EDIT 08/11/2014
% Thresholding of background added, in case user does not specify
%
% EDIT 08/13/2014
% Convergence criteria modified to be relative change in cost per voxel (to avoid weird behavior if cost goes below zero
% Flexible number of gaussians for equidistant means and same variance (N3 like) now works
% Function header added
%
% EDIT 08/20/2014
% (preliminary) support for priors added
%
% EDIT 02/10/2014
% Additional testing performed - release v. 0.9.0
%
% EDIT 10/13/2014
% various cleanups - release v. 0.9.1
%
% EDIT 10/14/2014
% correction using priors now work
% additional figures shown of priors and posteriors (implicit segmentation)
%
% EDIT 10/15/2014
% Now updated to support Matlab 2014b
%
% EDIT 10/17/2014
% Priors further tested, some cleanups and re-partitioning of code into separate functions header text updated - release v. 0.9.2
%
% EDIT 10/24/2014
% Some cleanups and bugfixes, a few functions made into subfunctions, preliminary MRF prior support - release v. 0.9.3
%
% EDIT 11/05/2014
% Corrected a minor bug in regularization matrix computation. Added a check for zero-intensity voxels in final mask. - release v. 0.9.4
%
% EDIT 11/12/2014
% Added supervoxel functionality, minor cleanups - version error check re-introduced
%
% EDIT 11/13/2014
% Re-added the N3 background thresholding, minor bug fixes - release v. 0.9.5
%
% EDIT 12/02/2014
% Great deal of changes: Removed most of the multi-resolution stuff for simplification (and it doesn't look like it's needed).
% Removed MRF mixture. Same variance mode greatly simplified. Various cleanups. - release v.0.9.6
%
% EDIT 04/09/2015
% Yet another big round of changes. Most notably supervoxels now work properly, and can also be used to filter away background.
% Longitudinal correction now also supported. - release v.0.9.9
%
% EDIT 08/15/2015
% Longitudinal bug wrt variance in difference field fixed (interdependency between bias field coefficients and difference variance).
% Many other fixes. - release v.1.0.0
%
% EDIT 09/05/2015
% Longitudinal bug fixed. - release v.1.0.0
%
% EDIT 12/13/2015
% Code MIT licensed. Public release (1.0.3)



if nargin < 1
    error('Call: [correctedData, OPTIONAL] = icc(data, OPTIONAL)')
end

versionString = strsplit(version, '.');
if (str2double([versionString{1}, '.', versionString{2}]) < 8.2)
    error('You need Matlab version 8.2 (R2013b) or higher for this program to run');
end


% =============================================================
%
% Parse inputs
%
% =============================================================
parser = inputParser;
addRequired(parser, 'data', @isnumeric);

addParameter(parser, 'mask', [], ...
    @(x) isnumeric(x) && all(size(x) == size(data)));
addParameter(parser, 'stepSize', [1,1,1], @isnumeric);
addParameter(parser, 'desiredStepSize', 1, @isnumeric);
addParameter(parser, 'backgroundThresholdType', 'none', ...
    @(x) any(validatestring(x, {'none', 'mixture', 'otzu', 'N3', 'supervoxel'})));
addParameter(parser, 'maxNumberIterations', 10000, ...
    @(x) isinteger(x) && x > 0);
addParameter(parser, 'maxNumberGMMIterations', 1000, ...
    @(x) isinteger(x) && x > 0);
addParameter(parser, 'prior', [], ...
    @(x) isnumeric(x) );
addParameter(parser, 'priorThreshhold', 0.99, ...
    @(x) isnumeric(x) );
addParameter(parser, 'superVoxelStepSize', 50, ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(parser, 'superVoxelTheta', 100, ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(parser, 'numberOfGaussians', 3, ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(parser, 'equidistantMeans', 0);
addParameter(parser, 'updateMeans', 1);
addParameter(parser, 'varianceInit', 0);
addParameter(parser, 'sameVariance', 0);
addParameter(parser, 'updateVariance', 1);
addParameter(parser, 'gmmFitType', 'mixture', ...
    @(x) any(validatestring(x, {'mixture', 'mixtureWithLabelPrior', 'mixtureWithSuperVoxel', 'wiener'})));
addParameter(parser, 'convergenceType', 'stdDifferenceField', ...
    @(x) any(validatestring(x, {'cost', 'stdDifferenceField'})));
addParameter(parser, 'convergence', 1e-5, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'convergenceGMM', 1e-6, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'filterPercentile', 0.00, ...
    @(x) isnumeric(x) && x > 0);


% todo: make 'bspline' option that uses proper spline scaling + regularization
% edit: not necessary, as long as regularization is scaled properly
% according to splines
addParameter(parser, 'smoothingType', 'bspline', ...
    @(x) any(validatestring(x, {'bspline', 'spm'})));
addParameter(parser, 'smoothingDistance', 50, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'smoothingRegularization', 1e-5, ...
    @(x) isnumeric(x));
addParameter(parser, 'voxelScaleRegularization', 1);
% addParameter(parser, 'smoothExponentiatedBias', 0);

addParameter(parser, 'fwhm', 0.15, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'wienerFilterNoise', 0.1, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'mimicN3', 'no',  ...
    @(x) any(validatestring(x, {'no', 'model', 'strict'})));

addParameter(parser, 'histogramBins', 200, ...
    @(x) isnumeric(x) && x > 0);
addParameter(parser, 'showPlots', 0);
addParameter(parser, 'verbosity', 'none', ...
    @(x) any(validatestring(x, {'none', 'info', 'debug'})));


parse(parser, data, varargin{:});

% data handling
mask = parser.Results.mask;
stepSize = parser.Results.stepSize;
desiredStepSize = parser.Results.desiredStepSize;
bkgThresType = parser.Results.backgroundThresholdType;

% model config
prior = parser.Results.prior;
priorThreshhold = parser.Results.priorThreshhold;
superVoxelStepSize = parser.Results.superVoxelStepSize;
% we define the theta as a 100%
superVoxelTheta = parser.Results.superVoxelTheta * 0.01;
numGaussians = parser.Results.numberOfGaussians;
equidistantMeans = parser.Results.equidistantMeans;
updateMeans = parser.Results.updateMeans;
sameVariance = parser.Results.sameVariance;
updateVariance =  parser.Results.updateVariance;
varianceInit = parser.Results.varianceInit;
gmmFitType = lower(parser.Results.gmmFitType);

% convergence settings
convergenceType = lower(parser.Results.convergenceType);
convergence = parser.Results.convergence;
convergenceGMM = parser.Results.convergenceGMM;
maxIters = parser.Results.maxNumberIterations;
maxMMParamIters = parser.Results.maxNumberGMMIterations;
filterPercentile = parser.Results.filterPercentile;

% bias shape
smoothingType = lower(parser.Results.smoothingType);
smoothingDist = parser.Results.smoothingDistance;
lambda = parser.Results.smoothingRegularization;
voxelScaleReg = parser.Results.voxelScaleRegularization;
% smoothExpBias = parser.Results.smoothExponentiatedBias;

% output
histogramBins = parser.Results.histogramBins;
showPlots = parser.Results.showPlots;
verbosity = lower(parser.Results.verbosity);

% N3 specific
fwhm = parser.Results.fwhm;
wienerFiltNoise = parser.Results.wienerFilterNoise;
mimicN3 = parser.Results.mimicN3;

if ~isempty(strfind(gmmFitType, 'mixturewithlabelprior')) || ~isempty(strfind(gmmFitType, 'mixturewithsupervoxel'))
    warning(['The model ', parser.Results.gmmFitType, ' is still highly experimental. It is supplied for testing purposes only.'])
end

% =================================================================
%
% Hard coded parameters - Settings we do not want to set outside the function right now
%
% =================================================================

% Moved to here, because it's really the only sensible setting, and the only time we want to deviate, is if doing strict N3
% interpolation when downsampling may mess up our priors, so we don't support it (and interpolation is bad for you in general)
downSamplingType = 'noInterp';

% for stepping through bias field smoothness using the regularizer
useMultiResolutionLambda = false;
lambdaStepValues = [9*lambda, 3*lambda, lambda];
currentLambda = lambda;
lambdaStep = 1;

convergenceSVMask = 1e-5;

minIntensityCap = 1;

%updateSpatialVariance = 1;
%updateWeights = 0;
updateSpatialVariance = 0;
updateWeights = 1;


% =============================================================
%
% Relevant checks and parameter dependent configs
%
% =============================================================

% determine if we're in longitudinal mode
T = 1;
numTP = size(data,4);
if numTP > 1
    A = eye(numTP);
    A(:,1) = 1/numTP;
    [Q,~] = qr(A);
    T = -Q';
end
%T = [0.5, 0.5; 1 -1];

% unless using a spatial prior, we assume one class consisting of multiple gaussians (or alternatively, several classes of one gaussian each)
numClasses = 1;

% if n3 flag is model or strict, ensure these options are set, since they cannot be changed in N3
if isempty(strfind(mimicN3, 'no'))
    if isempty(strfind(gmmFitType, 'wiener'))
        fprintf('Changing gaussian mixture fit mode to wiener\n');
        gmmFitType = 'wiener';
    end
    if isempty(strfind(smoothingType, 'bspline'))
        fprintf('Changing residual smoothing to N3 (bsplines)\n');
        smoothingType = 'bspline';
    end
    if numGaussians ~= 200
        fprintf('Changing numGaussians to 200\n');
        numGaussians = 200;
    end
    
    if isempty(strfind(convergenceType, 'stddifferencefield'))
        fprintf('Changing conversion type to standard deviation of difference (bias) field\n');
        convergenceType = 'stddifferencefield';
    end
    
    if ~voxelScaleReg
        fprintf('Changing scaling of regularization by number of voxels to true\n');
        voxelScaleReg = true;
    end
    
    % used in select cases
    n3Variance = 1 / (8 * log(2)) * fwhm^2;
    numLambdaSteps = 1;
end

% further set strict options for mimicN3
if ~isempty(strfind(mimicN3, 'strict'))
    % this stuff is simply too wacko, and therefore written out.
%     if ~smoothExpBias
%         fprintf('Changing smoothing of exponentiated bias estimate to true\n');
%         smoothExpBias = true;
%     end
    
    if isempty(strfind(downSamplingType, 'nearestneighbor'))
        fprintf('Changing downsampling type to nearestNeighbor\n');
        downSamplingType = 'nearestneighbor';
    end
    
    if isempty(mask) && isempty(strfind(bkgThresType, 'N3'));
        fprintf('Changing background thresholding type to N3\n');
        bkgThresType = 'N3';
    end
    
    % this is needed to 'translate' the regularization parameter to the 'real' regularization used by N3,
    % which is significantly higher than what you think when specifying it - because N3 does not multiply
    % by 2*variance properly - but we do, according to the derivations.
    lambda = lambda/(2*n3Variance);
    currentLambda = lambda;
end

% if the user didn't provide a spatial prior, masks or configured the background thresholding, set to default supervoxel
if isempty(prior) && isempty(mask) && ~isempty(strfind(bkgThresType, 'none'))
    if isempty(strfind(verbosity, 'none'))
        fprintf('Setting background Thresholding type to supervoxel.\n');
    end
    bkgThresType = 'supervoxel';
end

% wiener filtering does have a few dependencies...
switch gmmFitType
    case 'wiener'
        sameVariance = true;
        equidistantMeans = true;
        useParzenWindowHist = 1;
        useInterpTable = 1;
        doPad = 1;
    case 'mixturewithlabelprior'
        if isempty(prior)
            error('No spatial prior supplied');
        end
        
        % The numclasses parameter defines the number of classes we want to fit + a garbage class
        % verify that priors are only supported for the real classes
        numClasses = length(numGaussians);
        numPriors = size(prior, 2);
        if numClasses ~= numPriors + 1
            error('The vector numGaussians should have length equal to the number of priors + 1 (the background)');
        end
        
        % We need a map of which gaussians belong to what class, since we support multiple gaussians per class
        priorMap = zeros(1,sum(numGaussians));
        priorClass = 1;
        count = 1;
        for i = 1:numClasses
            for j = 1:numGaussians(i)
                priorMap(count) = priorClass;
                count = count + 1;
            end
            priorClass = priorClass + 1;
        end
        numGaussians = sum(numGaussians);
        
        % compute background prior + ensure they all sum up to one
        totalProbability = sum(prior, 2);
        idx = find(totalProbability > 1);
        entries = length(idx);
        if(entries)
            if isempty(strfind(verbosity, 'none'))
                fprintf('Probabilities sum to over one, probably a resampling artifact. Normalizing %d entries to ensure background probability does not become negative\n', entries);
            end
            for i = 1:numPriors
                prior(idx, i) = prior(idx, i) ./ totalProbability(idx);
            end
            totalProbability(idx) = 1;
        end
        prior(:,numPriors+1) = 1 - totalProbability;
    case 'mixturewithsupervoxel'
        bkgThresType = 'supervoxel';
end


% =============================================================
%
% Plotting stuff
%
% =============================================================

% figure configs
estimatedBiasFieldFig = -1;
correctedDataFig = -1;
historyFig = -1;
gmmFig = -1;
gmmCostFig = -1;

% for plotting stuff
if showPlots
    gmmFig = figure('name', 'gmm');
    gmmCostFig = figure('name', 'gmm cost');
    maskFig = figure('name', 'masks');
    estimatedBiasFieldFig = figure( 'name', 'bias' );
    correctedDataFig = figure( 'name', 'corrected' );
    historyFig = figure( 'name', 'History' );
    bFig = figure('name', 'basis functions and regularizer');
    if ~isempty(prior)
        priorFig = figure('name', 'priors');
        posteriorFig = figure('name', 'posteriors');
    end
    if numTP
        diffVarFig = figure('name', 'difference variance');
        diffHistFig = figure('name', 'difference Histogram');
    end
end

% =================================================================
%
% Preprocessing
%
% =================================================================

% Time it
timeFitData = 0;
timeEstBias = 0;
timeFitWeights = 0;
startTime = tic;
startTimePreProcessing = tic;


% =================================================================
%
% Downsampling, zero intensity masking + log transform
%
% =================================================================

% save original data for final correction
origData = data;
origDataSize = size(origData);
origDataSize = origDataSize(1:3);
origStepSize = stepSize;
numTotalOrigVoxels = numel(origData(:,:,:,1));


% downsample according to desired resolution and sampling type
[data1, start, stepSize, sampleIDX] = downSampleData(origData(:,:,:,1), origStepSize, desiredStepSize, downSamplingType);
dataSize = size(data1);
numTotalVoxels = numel(data1);
data = zeros([dataSize, numTP]);
data(:,:,:,1) = data1;
clear data1;

% filter away top percentile of voxels (garbage/artifacts)
filterBins = max(data(:));
maxIntensityCap = filterBins + eps;
[histogram, binCenters] = hist(data(:), ceil(filterBins));
histogram = histogram / sum(histogram(:));
cumHist = cumsum(histogram);
[~, idx] = find(cumHist > (1-filterPercentile));
if ~isempty(idx)
    maxIntensityCap = binCenters(idx(1));
end

% non zero masks (needed for log transformation)
minMask = data(:,:,:,1) > minIntensityCap;
maxMask = data(:,:,:,1) < maxIntensityCap;

% note that we do not include a maxOrigMask - this is intended, to avoid left out voxels in e.g., white matter after correction
minOrigMask = origData(:,:,:,1) > minIntensityCap;
%maxOrigMask = origData(:,:,:,1) < maxIntensityCap;

% longitudinal
for i = 2:numTP
    [data(:,:,:,i)] = downSampleData(origData(:,:,:,i), origStepSize, desiredStepSize, downSamplingType);  
    minOrigMask = minOrigMask .* (origData(:,:,:,i) > minIntensityCap);
  %  maxOrigMask = maxOrigMask .* (origData(:,:,:,i) < maxIntensityCap);
    minMask = minMask .* (data(:,:,:,i) > minIntensityCap);
    maxMask = maxMask .* (data(:,:,:,i) < maxIntensityCap);
end

% logical masks
minMaskLogical = minMask > 0;
minMaxMaskLogical = minMask > 0 & maxMask > 0;
minOrigMaskLogical = minOrigMask > 0;

% need these for later
dataBackup = data;
% using the original data, rather than zeros, prevents from problems due to improper bi-modal masking
% e.g., contrast bias results in missing temporal lopes
% this way we preserve image information, even if the intensities were not corrected due to the masking
correctedData = origData; 

% log transform the data
for i = 1:numTP  
    tmpData = data(:,:,:,i);
    tmpData(minMaskLogical) = log(tmpData(minMaskLogical));
    data(:,:,:,i) = tmpData .* minMask;
    
    tmpOrigData = origData(:,:,:,i);
    tmpOrigData(minOrigMaskLogical) = log(tmpOrigData(minOrigMaskLogical));
    origData(:,:,:,i) = tmpOrigData .* minOrigMask;
end
clear tmpData tmpOrigData

% Transform it to feature space (1 dimension is the common signal, remaining ones are noise + difference bias
data = reshape((T * reshape(data, [numTotalVoxels, numTP])')', [dataSize, numTP]);
origData = reshape((T * reshape(origData, [numTotalOrigVoxels, numTP])')', [origDataSize, numTP]);

% signal is the real biology (response from tissue)
signalData = data(:,:,:,1);

% trust criterion used for gaussians
% Trust criterions that ensures gaussians do not collapse. A value of 1e-2 corresponds to 1% for e.g., the posterior contribution
trustWeights = 1e-3;%trustCriterion / numGaussians;
%this is a bit ad hoc...
trustVariance = (0.01*(max(signalData(minMaxMaskLogical))-min(signalData(minMaxMaskLogical))))^2;


% =================================================================
%
% Initialize parameters for super voxel mixture (or superVoxel masking)
%
% =================================================================

if ~isempty(strfind(gmmFitType, 'mixturewithsupervoxel')) || ~isempty(strfind(bkgThresType, 'supervoxel'))
    % The number of supervoxels along each dimension is given by the superVoxel spacing parameter in mm
    numSuperVoxels = round(origDataSize .* origStepSize / superVoxelStepSize);
    
    % we need to store the positions of all voxels
    positionsX = repmat(repmat(((1:dataSize(1))' - 0.5), [1, dataSize(2)]), [1, 1, dataSize(3)]);
    positionsY = repmat(repmat(((1:dataSize(2)) - 0.5), [dataSize(1), 1]), [1, 1, dataSize(3)]);
    positionsZ = reshape(ones(dataSize(1)*dataSize(2), 1) * ((1:dataSize(3)) - 0.5), [dataSize(1), dataSize(2), dataSize(3)]);
    
    dataPositions = [positionsX(:), positionsY(:), positionsZ(:)];
    dataPositions = dataPositions .* repmat(stepSize, [numTotalVoxels, 1]);
 
    % initialize supervoxel centers
    svWidthX = dataSize(1) / numSuperVoxels(1);
    svPosX = round(svWidthX * (0.5 + (0:(numSuperVoxels(1)-1))));
    svPosX = repmat(repmat(svPosX', [1, numSuperVoxels(2)]), [1, 1, numSuperVoxels(3)]);
    svPosX = svPosX(:);
    
    svWidthY = dataSize(2)  / numSuperVoxels(2);
    svPosY = round(svWidthY * (0.5 + (0:(numSuperVoxels(2)-1))));
    svPosY = repmat(repmat(svPosY, [numSuperVoxels(1), 1]), [1,1,numSuperVoxels(3)]);
    svPosY = svPosY(:);
    
    svWidthZ = dataSize(3) / numSuperVoxels(3);
    svPosZ = round(svWidthZ * (0.5 + (0:(numSuperVoxels(3)-1))));
    svPosZ = reshape(ones(numSuperVoxels(1)*numSuperVoxels(2), 1) * svPosZ, [numSuperVoxels(1), numSuperVoxels(2), numSuperVoxels(3)]);
    svPosZ = svPosZ(:);
    
    superVoxelCenters = [svPosX(:), svPosY(:), svPosZ(:)];
    clear svPosX svPosY svPosZ positionsX positionsY positionsZ

    % One gaussian per superVoxel
    numGaussiansSV = prod(numSuperVoxels);
end


% =================================================================
%
% Masking
%
% =================================================================

% if no mask or prior has been specified, we need to mask out the background
bkgThres = 0;
origMaskBKG = ones(origDataSize);
maskBKG = ones(dataSize);
if isempty(strfind(bkgThresType, 'none'))
    if isempty(strfind(verbosity, 'none'))
        fprintf('Creating mask by tresholding background\n');
    end
   
    % I've removed the orig data thresholding, since it makes no sense to do it twice; we need only ONE threshold, and we want to be
    % sure that the same voxels, downsampled or not, are either in both masks - or one. Computing twice potentially messes this up
    switch bkgThresType
        case {'otzu', 'n3', 'mixture'}
            for i = 1:numTP
                % the point of this is to work on the non-log, non-feature transformed data (if we so want to), 
                % whereas the supervoxel stuff requires log 
                tmpData = dataBackup(:,:,:,i);
                % correctedData is just a backup for now
                tmpOrigData = correctedData(:,:,:,i);
                
                % oook. Use a histogram that bins intensities using the range to define the number of bins
                nBins = ceil(max(tmpData(:))-min(tmpData(:)) + 1);

                % note that N3 actually does some sort of rounding on the subsampled threshing, but not in the end when threshing
                % the full res data. I was unable to reproduce the exact value they arrive at, so this is a (very) minor deviation
                % from how N3 does things
                % BTW, N3 does the threshing on the downsampled data, and AGAIN at full res in the end. I've removed this...
                bkgThres = biModalThreshold(tmpData, false, bkgThresType, nBins);

                maskBKG = maskBKG .* (tmpData >= bkgThres);
                origMaskBKG = origMaskBKG .* (tmpOrigData >= bkgThres);  
            end         
        case 'supervoxel'    
            origMaskBKG = zeros(origDataSize);
            % pull out the common signal
            
            tmpData = data(:,:,:,1);
            tmpOrigData = origData(:,:,:,1);
            tmpData = [tmpData(minMaxMaskLogical), dataPositions(minMaxMaskLogical,:)];
           % tmpOrigData = tmpOrigData(nonZeroOrigMaskLogical);
                 
            % some values that we will need for initialization of the mixture fit later
            means = repmat(mean(tmpData(:,1)), [numGaussiansSV, 1]);
            means = [means, superVoxelCenters];
            
            % default variance initialization
            variances = repmat(var(tmpData(:,1)), [numGaussiansSV, 1]);
            variances = [variances, repmat(superVoxelStepSize^2, [numGaussiansSV, 3])]; 
            
            % just use flat weights...
            weights = ones(numGaussiansSV, 1) / numGaussiansSV;
            
            % do one full round of supervoxel fitting to the image, in order to get a supervoxel segmentation
            % that can be used to do a background/foreground separation
            [means, variances, weights] = fitGMMwithSV(...
                tmpData, ...
                means, ...
                variances, ...
                weights, ...
                superVoxelTheta, ...
                1, ...
                updateSpatialVariance, ...
                updateWeights, ...
                0, ...
                maxMMParamIters, ...
                convergenceSVMask, ...
                trustWeights, ...
                trustVariance, ...
                histogramBins, ...
                verbosity, ...
                gmmFig, ...
                gmmCostFig);
            
            % run again, but just once to obtain posterior for all voxels
            % do one full round of supervoxel fitting to the image, in order to get a supervoxel segmentation
            % that can be used to do a background/foreground separation we need to store the positions of all voxels
            % note that we do not include the maxOrigMask - this is intended, to avoid left out voxels in e.g., white matter
            
            positionsX = repmat(repmat(((1:origDataSize(1))' - 0.5), [1, origDataSize(2)]), [1, 1, origDataSize(3)]);
            positionsY = repmat(repmat(((1:origDataSize(2)) - 0.5), [origDataSize(1), 1]), [1, 1, origDataSize(3)]);
            positionsZ = reshape(ones(origDataSize(1)*origDataSize(2), 1) * ((1:origDataSize(3)) - 0.5), [origDataSize(1), origDataSize(2), origDataSize(3)]);
            
            origDataPositions = [positionsX(:), positionsY(:), positionsZ(:)];
            clear positionsX positionsY positionsZ 
            
            % one gmm iteration to obtain the posteriors for all voxels
            origDataPositions = origDataPositions .* repmat(origStepSize, [numTotalOrigVoxels, 1]);
            
            voxelSegMeans = [];
            numelPerSlize = origDataSize(1)*origDataSize(2);
            for i = 1:origDataSize(3)

                tmpSliceData = tmpOrigData(:,:,i);
                sliceData = [tmpSliceData(:), origDataPositions((1+(i-1)*numelPerSlize):(i*numelPerSlize),:)];
                sliceMask = minOrigMaskLogical(:,:,i);

                [~, ~, ~, ~, posteriors] = fitGMMwithSV(...
                    sliceData(sliceMask, :), ...
                    means, ...
                    variances, ...
                    weights, ...
                    superVoxelTheta, ...
                    1, ...
                    updateSpatialVariance, ...
                    updateWeights, ...
                    0, ...
                    1);

                [~,voxelSVSeg] = max(posteriors,[], 2);
                voxelSegMeans = [voxelSegMeans; means(voxelSVSeg, 1)]; %#ok;
            end
          
            clear posteriors voxelSVSeg sliceData sliceDataPositions sliceMask numelPerSlize

            % do a bimodal filtering of supervoxel center intensity means, given the previous fit
            bkgThres = biModalThreshold(means(:,1), 0, 'otzu', numGaussiansSV);
            if ~isempty(strfind(verbosity, 'debug'))
                fprintf('background Threshold %d\n', bkgThres);
            end  
            idx = find(means(:, 1) > bkgThres);   
            
            % serves as initialization for later - makes first iteration super fast
            means = means(idx, :);
            variances = variances(idx, :);
            weights = weights(idx);
            % necessary to maintain proper mixture
            weights = weights / sum(weights);

            % compute a mask, given the hard seg'ed voxels
            origMaskBKG(minOrigMaskLogical) = voxelSegMeans(:) > bkgThres;       
            [maskBKG, ~, ~, ~] = downSampleData(origMaskBKG, origStepSize, desiredStepSize, downSamplingType);
    end
end
clear tmpData tmpOrigData

% if provided, downsample mask
origMask = ones(origDataSize);
if isempty(mask)
    mask = ones(dataSize);
else
    for i = 1:numTP
        origMask = origMask .* mask(:,:,:,i);
    end
    mask = downSampleData(origMask, origStepSize, desiredStepSize, downSamplingType);
end

% if prior provided, downsample prior mask
if isempty(prior) || numel(prior) == 1
    origMaskPrior = ones(origDataSize);
    maskPrior = ones(dataSize);
else
    origMaskPrior = zeros(origDataSize);
    % mask out voxels which are threshhold % or more likely to be garbage/background + voxels that are virtually zero intensities
    origMaskPrior(prior(:, end) < priorThreshhold  & origData(:) > minIntensityCap) = 1;
    
    % ok as long as we only allow downsampling where values are not interpolated
    maskPrior = downSampleData(origMaskPrior, origStepSize, desiredStepSize, downSamplingType);
end

% now we compute the full mask
% note that we do not include the maxOrigMask - this is intended, to avoid left out voxels in e.g., white matter
origMask = origMask .* origMaskPrior .* origMaskBKG .* minOrigMask;
mask = mask .* maskPrior .* maskBKG .* minMask .* maxMask;

if showPlots
    figure(maskFig);
    subplot(2,2,1); showImage(mask);
    subplot(2,2,2); showImage(minMask);
    subplot(2,2,3); showImage(maskBKG);
    subplot(2,2,4); showImage(maskPrior);
end

% we use logical indexing heavily, so this saves us from quite a lot of repeated boolean operations
origMaskLogical = origMask > 0.5;
maskLogical = mask > 0.5;

% if we're using priors, we need to filter out all voxels that are not used, after downsampling
if ~isempty(prior) && numel(prior) > 1
    priorsFiltered = prior(sampleIDX, :);
    priorsFiltered = priorsFiltered(maskLogical, :);
    
    if showPlots
        figure(priorFig)
        numFigsPerDim = ceil(sqrt(numClasses));
        for i  = 1:numClasses
            subplot(numFigsPerDim,numFigsPerDim,i);
            tmp = mask;
            tmp(mask > 0.5) = priorsFiltered(:,i);
            showImage(tmp);
        end
    end
end

% This way we know exactly how many voxels we operate on when computing the minLogLikelihood
numVoxels = sum(maskLogical(:));

if(numVoxels ~= sum(mask(:)))
    error('Mask has not been properly binarized');
end

% Only work on non-zero intensities
intensities = signalData(maskLogical);

% =================================================================
%
% Initialize parameters for one class mixture and priors
%
% =================================================================

% if we're running the supervoxel model, these have already been initialized when creating the mask
switch gmmFitType
    case 'mixturewithsupervoxel'
        dataPositions = dataPositions(maskLogical, :);
    case 'mixturewithlabelprior'
        means = zeros(1,numGaussians);
        variances = zeros(1,numGaussians);
        weights = zeros(1,numGaussians);
        
        % compute class statistics
        posteriors = priorsFiltered;
        posteriorsSum = sum(posteriors);
        
        % update weights - flat initially for gaussians within each class
        for l = 1:numClasses
            weights(priorMap == l) = 1 / sum(priorMap == l);
        end
        
        % update class means and variances
        classMeans = sum(posteriors .* repmat(intensities,[1, numClasses])) ./ (posteriorsSum + eps);
        squaredDists = zeros(numVoxels, numClasses);
        for i = 1:numClasses
            squaredDists(:, i) = (intensities - classMeans(i)).^2;
        end
        classVariances = sum(posteriors .* squaredDists) ./ (posteriorsSum + eps);
        classStd = sqrt(classVariances);
        
        % now that we have some per class statistics, we can initialise the gaussians within each
        for i = 1:numClasses
            gaussiansInClass = sum(priorMap == i);
            range = [classMeans(i)-3*classStd(i), classMeans(i)+3*classStd(i)];
            rangeSpan = range(2) - range(1);
            gaussianWidth = rangeSpan / (gaussiansInClass - 1);
            
            means(priorMap == i) = range(1) + (0:gaussiansInClass-1) * gaussianWidth;
            
            variance = gaussianWidth^2;
            if varianceInit
                variance = varianceInit;
            end
            variances(priorMap == i) =  repmat(variance, [1, gaussiansInClass]);
        end
            % one class configurations (mixture, N3)
    case {'mixture','wiener'}     
        range = [min(intensities(:)), max(intensities(:))];
        % Define range of intensities that our non-zero data spans (same as DHistogram.cc in N3)
        rangeSpan = range(2) - range(1);
        gaussianWidth = rangeSpan / (numGaussians - 1);
        
        % mean init
        means = range(1) + (0:numGaussians-1)*gaussianWidth;
        
        % default variance initialization
        variance = gaussianWidth^2;
        if varianceInit
            variance = varianceInit;
        end
        variances = repmat(variance, [1, numGaussians]);
        
        
        % uniform weights - tests indicate that using the histogram for init doesn't really make any difference
        % well, apparently this is not always the case... in the same variance scenario, using the histogram seems to be pretty decent.
        % If I do not do it, the method becomes more sensitive to how I scale the variance - which MIGHT perform comparably, but not always
        weights = ones(1, numGaussians) / numGaussians;
        if sameVariance
            weights = hist(intensities, numGaussians) / sum(hist(intensities, numGaussians));
        end
end


% =================================================================
%
% Initialize histories
%
% =================================================================

oneOverEps = 1 / eps;

costHist = oneOverEps;
dataTermHist = oneOverEps;
regTermHist = oneOverEps;

stdFieldHist  = oneOverEps;
stdDiffHist  = oneOverEps;
cvDataHist  = oneOverEps;

meansHist = means;
variancesHist = variances;
weightsHist = weights;
infoStruct.bkgThresh = bkgThres;

% if useMultiResolutionLambda
%     lambdaHist = 1;
% end

% =================================================================
%
% Initialize parameters for bias field
%
% =================================================================


% pre-allocate
estimatedBiasField = zeros(dataSize);
weightImage = zeros(dataSize);
residualImage = zeros(dataSize);

% compute everything needed for basis function bias field representation
[numBFuncs, domain] = computeNumberOfSeparableBasisFunctions(origDataSize, origStepSize, smoothingType, smoothingDist);
[B, B2, indexMapB2] = computeSeparableBasisFunctions(numBFuncs, domain, smoothingDist, dataSize, smoothingType,  stepSize, start);
J = computeRegularizers(smoothingType, numBFuncs, domain, origStepSize);

% speed-up, A2 does not depend on variance
if sameVariance
    A2 = constructAtWA(mask, B2, indexMapB2);
end

% plot the basis functions + regularizer
if showPlots
    figure(bFig);
    subplot(2,2,1); plot(B{1});     axis tight;
    subplot(2,2,2); plot(B{2});     axis tight;
    subplot(2,2,3); plot(B{3});     axis tight;
    subplot(2,2,4); imagesc(J);     axis tight;
end

% initialize coefficients (flat field in log-domain)
C = zeros(prod(numBFuncs), numTP);

% Balance the data term and the reg term (corresponds to dividing the data term by N)
if voxelScaleReg
    J = J * numVoxels;
end

convDiffField = signalData(maskLogical);

% end of preprocessing
timePreProcessing = toc(startTimePreProcessing);


% print some info
if isempty(strfind(verbosity, 'none'))
    fprintf('Settings:\n');
    fprintf('==================================\n');
    fprintf('Number of time points: %d\n', numTP);    
    fprintf('Number of Gaussians: %d\n', numGaussians);
    fprintf('Equidistant Gaussian mean restriction: %d\n', equidistantMeans);
    fprintf('Same variance Gaussian restriction: %d\n', sameVariance);
    fprintf('Update Gaussian variance during estimation: %d\n', updateVariance);
    fprintf('Gaussian variance initialization scale: %.3f\n', varianceInit);
    fprintf(['Type of Mixture model fitting: ', gmmFitType, '\n']);
    fprintf(['Convergence scheme: ', convergenceType, '\n']);
    fprintf('Stop criterion: %g\n', convergence);
    fprintf(['Smoothing type: ', smoothingType, '\n']);
    fprintf('Smoothing distance: %d\n', smoothingDist);
    fprintf('lambda: %g\n', lambda);
    fprintf('Scale regularization by number of voxels: %d\n', voxelScaleReg);
    fprintf('Histogram bins: %d\n', histogramBins);
    fprintf('==================================\n\n');
end

if ~isempty(strfind(verbosity, 'debug'))
    fprintf('Bias field correcting data:\n\n');
end

% =================================================================
%
% Bias field correct data
%
% =================================================================

% iterative EM scheme starts here, iteration cap is given by maxIters
for iter = 1:maxIters
    if ~isempty(strfind(verbosity, 'debug'))
        fprintf('ITERATION %d\n', iter);
    end
    
    
    % bias field is always estimated with respect to original data; only in the end is the data corrected.
    biasFieldCorrectedData = signalData - estimatedBiasField;
    convDiffField = convDiffField - estimatedBiasField(maskLogical);
    
    if showPlots
        figure(correctedDataFig);
        showImage( exp(biasFieldCorrectedData) .* mask);
        figure(estimatedBiasFieldFig);
        showImage(exp(estimatedBiasField) .* mask);
    end
    
    % N3 doesn't even use the coefficient of variation for termination - it uses the standard deviation!
    % but ok, the CV is just a normalization...
    stdFieldHist(end + 1) = std(estimatedBiasField(maskLogical)); %#ok
    stdDiffHist(end + 1) = std(convDiffField); %#ok
    % coefficient of variation within mask on subsampled data
    cvDataHist(end + 1) = std(biasFieldCorrectedData(maskLogical)) / mean(biasFieldCorrectedData(maskLogical)); %#ok
    
    % pull out the current set of bias field corrected intensities for next mixture model fit
    intensities = biasFieldCorrectedData(maskLogical);
    
   % dataTermHist(end) + currentLambda * C(:,1)' * J * C(:,1)
    
    startTimeFitData = tic;
    switch gmmFitType
        % N3 fit
        case 'wiener'
            [means, variances, weights, minLogL, expectation, elapsedTimeWeightsIter] = fitGMMWiener(...
                intensities, ...
                numGaussians, ...
                fwhm, ...
                wienerFiltNoise,...
                histogramBins, ...
                useParzenWindowHist, ...
                useInterpTable, ...
                doPad, ...
                gmmFig);
            timeFitWeights = timeFitWeights + elapsedTimeWeightsIter;
            
            % regular GMM fit
        case 'mixture'
            [means, variances, weights, minLogL, ~, Si, expectation] = fitGMM(...
                intensities, ...
                means, ...
                variances, ...
                weights, ...
                equidistantMeans, ...
                updateMeans, ...
                sameVariance, ...
                updateVariance, ...
                maxMMParamIters, ...
                convergenceGMM, ...
                trustWeights, ...
                trustVariance, ...
                histogramBins, ...
                verbosity, ...
                gmmFig, ...
                gmmCostFig);
            % regular GMM fit + priors
        case 'mixturewithlabelprior'
            [means, variances, weights, minLogL, posteriors, Si, expectation] = fitGMMwithPrior(...
                intensities, ...
                priorsFiltered, ...
                priorMap, ...
                means, ...
                variances, ...
                weights, ...
                maxMMParamIters, ...
                convergenceGMM, ...
                trustWeights, ...
                trustVariance, ...
                histogramBins, ...
                verbosity, ...
                gmmFig, ...
                gmmCostFig);
            
            if showPlots && numel(prior) > 1
                figure(posteriorFig);
                segmentation = zeros(dataSize);
                for i = 1:numClasses
                    segmentation(maskLogical) = sum(posteriors(:, priorMap == i), 2);
                    subplot(2,2,i)
                    showImage(segmentation);
                end
            end
        case  'mixturewithsupervoxel'
            [means, variances, weights, minLogL, ~, Si, expectation] = fitGMMwithSV(...
                [intensities, dataPositions], ...
                means, ...
                variances, ...
                weights, ...
                superVoxelTheta, ...
                updateVariance, ...
                updateSpatialVariance, ...
                updateWeights, ...
                sameVariance, ...
                maxMMParamIters, ...
                convergenceGMM, ...
                trustWeights, ...
                trustVariance, ...
                histogramBins, ...
                verbosity, ...
                gmmFig, ...
                gmmCostFig);
        otherwise
            error('invalid Gaussian Mixture Model Fitting');
    end
    
    % Time the histogram fit
    timeFitData = timeFitData + toc(startTimeFitData);
    
    meansHist(:,:,end+1) = means; %#ok
    variancesHist(:,:,end+1) = variances; %#ok
    weightsHist(:,:,end+1) = weights; %#ok

    dataTermHist(end + 1) = minLogL; %#ok
    regTermHist(end + 1) = currentLambda * C(:,1)' * J * C(:,1);%#ok
    costHist(end + 1) = dataTermHist(end) + regTermHist(end);%#ok

    
    if ~isempty(strfind(verbosity, 'debug'))
        fprintf('\tCost: %.5f\n', costHist(end));
    end
    
    if showPlots
        figure(historyFig);
        % subplot(2,2,1);
        plot(costHist(2:end), 'b');
        axis tight
     %   hold on
       %   plot(dataTermHist(2:end), 'r');
      %    plot(regTermHist(2:end), 'g');
     %   hold off
        title('Cost');
        xlabel('Iteration');
        ylabel('Cost');
    end
    
    switch convergenceType
        case 'cost'
            relativeChangeCost = (costHist(end-1) - costHist(end)) / numVoxels;
        case 'stddifferencefield'
            relativeChangeCost = stdDiffHist(end);
        otherwise
            error('Convergence type not specified');
    end

    % Convergence / multi resolution
    if (abs(relativeChangeCost) < (convergence * 10))
        % Start with a very rigid field, and increase flexibility as we go along.
        % Avoid local minima when the field is too flexible compared to how well we're able to fit the data
        if useMultiResolutionLambda && (lambdaStep < numLambdaSteps)
            % increase smoothing flexibility by decreasing lambda (regularization)
            lambdaStep = lambdaStep + 1;
            currentLambda = lambdaStepValues(lambdaStep);
            lambdaHist(end + 1) = iter; %#ok
            
            if ~isempty(strfind(verbosity, 'debug'))
                fprintf('Decreasing lambda to %f \n', currentLambda);
            end
            % we're done!
        elseif (abs(relativeChangeCost) < convergence)
            if ~isempty(strfind(verbosity, 'debug'))
                fprintf('Converged in iteration %d!\n', iter);
            end
            break;
        end
    end
    
    % store previous fild for difference field computation in next iter
    convDiffField = estimatedBiasField(maskLogical);
    
    % Estimate and time the bias field
    startTimeEstBias = tic;
    
    % the noisy residual that we want smoothed
    residualImage(maskLogical) = signalData(maskLogical) - expectation;
    
    % fit basis functions to data to obtain coefficients
    if(sameVariance)
        lhs = A2 + variances(1) * 2 * currentLambda * J;
        rhs = constructAtWr(residualImage, mask, B);
    else
        weightImage(maskLogical) = Si;
        
        lhs = constructAtWA(weightImage, B2, indexMapB2) + 2 * currentLambda * J;
        rhs = constructAtWr(residualImage, weightImage, B);
    end
    C(:,1) = lhs \ rhs;
    
    % evaluate basis functions
    estimatedBiasField = expandSeparableBasisFunctions(C(:,1), B);
    
    %     if ~isempty(strfind(verbosity, 'debug'))
    %         fprintf('\tResidual min: %f, max: %f\n',  exp(min(residualImage(maskLogical))), exp(max(residualImage(maskLogical))));
    %         fprintf('\tField    min: %f, max: %f\n',  exp(min(estimatedBiasField(maskLogical))), exp(max(estimatedBiasField(maskLogical))));
    %     end
    timeEstBias = timeEstBias + toc(startTimeEstBias);
end

startTimeFitData = tic;

% longitudinal correction - may not be necessary, but here for simplicity
differenceVariance = zeros(numTP-1,1);
A2 = constructAtWA(mask, B2, indexMapB2);

% we're using regulariation, which makes difference variance and
% coefficients interdependent, so things are not longer closed form :-/
if lambda == 0
    maxIters = 1;
end

for i = 2:numTP
    tmp = data(:,:,:,i);
    figure
    if showPlots
        [histogram, binCenters] = hist(tmp(maskLogical), histogramBins);
        histogram = histogram / sum(histogram(:));
        binWidth = binCenters(2) - binCenters(1);        
        gaussPoints = linspace(min(tmp(maskLogical)), max(tmp(maskLogical)),512);
    end
    
    diffCostHist = 1/eps;
    biasV = zeros(dataSize); 
    for difIter = 1:maxIters
        % update variance directly (E-step)
        varianceV = sum((tmp(maskLogical)- biasV(maskLogical)).^2)/numVoxels;
        
        % convergence
        posteriors = (1 / sqrt(2*pi*varianceV)) * exp(-0.5 * tmp(maskLogical).^2 / varianceV);
        minLogL = -sum(log(sum(posteriors) + eps));
        diffCostHist(end + 1) = minLogL;%#ok      
        if (abs(diffCostHist(end-1) - diffCostHist(end)))/numVoxels < 1e-15
            if ~isempty(strfind(verbosity, 'debug'))
                fprintf('Difference estimation converged in iteration %d for image %d!\n', difIter, i);
            end
            break;
        end
        
        lhs = A2 + 2 * varianceV * lambda * J;
        rhs = constructAtWr(tmp, mask, B);
        C(:,i) = lhs \ rhs;
        biasV = expandSeparableBasisFunctions(C(:,i), B);
        
        if showPlots && ishandle(diffVarFig)
            figure(diffHistFig)
            gaussianEvals =  (1 / sqrt(2*pi*varianceV)) * exp(-0.5 * gaussPoints.^2 / varianceV) * binWidth;         
            clf
            grid
            bar(binCenters, histogram);
            hold on
            plot(gaussPoints, gaussianEvals, 'g', 'lineWidth', 1);
            xlabel('Bins')
            xlim([1,7]);
            title('Log-intensity difference histogram');
            ylabel('Normalized Voxel Count');
            set(gca, 'fontsize', 6)
            axis tight
            hold off

            figure(diffVarFig)
            plot(diffCostHist(2:end), 'b');
            axis tight
            title('Difference Variance Cost');
            xlabel('Iteration');
            ylabel('Variance');
            drawnow
        end
    end
    differenceVariance(i-1) = varianceV;
end
  
% Time the histogram fit
timeFitData = timeFitData + toc(startTimeFitData);

% Reconstruct the bias field at the original resolution
if ~isempty(strfind(verbosity, 'debug'))
    disp('Reconstructing the bias field at the original resolution')
end

startTimePostProcessing = tic;

% Transform Back
C = (T \ C')';


% % wacko N3 smoothing of bias in original domain.
% if smoothExpBias
%     residualImage(maskLogical) = exp(commonBiasField(maskLogical));
%     
%     % fit basis functions to data to obtain coefficients
%     if(sameVariance)
%         lhs = A2 + variances(1) * 2 * lambda * J;
%         rhs = constructAtWr(residualImage, mask, B);
%     else
%         lhs = constructAtWA(weights, B2, indexMapB2) + 2 * lambda * J;
%         rhs = constructAtWr(residualImage, weightImage, B);
%     end
%     
%     C(:,1) = lhs \ rhs;
% end

% compute basis functions at original resolution
B = computeSeparableBasisFunctions(numBFuncs, domain, smoothingDist, origDataSize, smoothingType, origStepSize, start);
estimatedBiasField = zeros([origDataSize, numTP]);
for i = 1:numTP
    % evaluate basis functions obtain bias at original resolution   
    tmpBiasField = expandSeparableBasisFunctions(C(:,i), B);
    % transform back
    tmpBiasField(origMaskLogical) = exp(tmpBiasField(origMaskLogical));
    % mask it
    estimatedBiasField(:,:,:,i) = tmpBiasField .* origMask;
    % recall correctedData contains the original Data  
    correctedData(origMaskLogical) = correctedData(origMaskLogical) ./ (tmpBiasField(origMaskLogical) + eps);
end

%for output
mask = origMask;

timePostProcessing = toc(startTimePostProcessing);

% Report on computation time
elapsedTime = toc(startTime);
oneOverElapsedTime100 = 100 / elapsedTime;
timePerIter = (timeFitData + timeEstBias) / (iter-1);

if isempty(strfind(verbosity, 'none'))
    disp('Done computing bias fields, displaying statistics');
    disp(['   Iterations used           : ', num2str(iter-1)]);
    disp(['   Time taken                : ', num2str(elapsedTime), ' seconds']);
    disp(['   Time per iteration        : ', num2str(timePerIter), ' seconds']);
    disp(['   Pre processing            : ', num2str(timePreProcessing * oneOverElapsedTime100) '%']);
    disp(['   Fitting models            : ', num2str(timeFitData * oneOverElapsedTime100) '%']);
    disp(['   Updating weights          : ', num2str(timeFitWeights * oneOverElapsedTime100) '%']);
    disp(['   Computing bias fields     : ', num2str(timeEstBias * oneOverElapsedTime100) '%']);
    disp(['   Post processing           : ', num2str(timePostProcessing * oneOverElapsedTime100) '%']);
end

if nargout > 3
    infoStruct.iterations = iter-1;
    infoStruct.elapsedTime = elapsedTime;
    infoStruct.timePerIteration = timePerIter;
    infoStruct.costHistory = costHist;
    infoStruct.dataTermCostHistory = dataTermHist;
    infoStruct.regularizerTermCostHist = regTermHist;
    infoStruct.cvDataHistory = cvDataHist;
    infoStruct.stdFieldHistory = stdFieldHist;
    infoStruct.stdDifferenceFieldHistory = stdDiffHist;
    infoStruct.meanHistory = meansHist;
    infoStruct.varianceHistory = variancesHist;
    infoStruct.differenceVariance = differenceVariance;
    infoStruct.weightHistory = weightsHist;
    infoStruct.preProcessing = timePreProcessing * oneOverElapsedTime100;
    infoStruct.fittingModels = timeFitData * oneOverElapsedTime100;
    infoStruct.computingBiasFields = timeEstBias * oneOverElapsedTime100;
    infoStruct.postProcessing = timePostProcessing * oneOverElapsedTime100;
    infoStruct.bkgThresh = bkgThres;
end

% write out segmentations
if nargout > 4
    segmented = zeros(origDataSize);
    switch gmmFitType
%         % regular GMM fit
%         case 'mixture'
% 
%          % regular GMM fit + priors
%         case 'mixtureWithLabelPrior'
%             % for mem efficiency...
%             numelPerSlize = origDataSize(1)*origDataSize(2);
%             for i = 1:origDataSize(3)       
%                 sliceData = origData(:,:,i);
%                 sliceMask = origMaskLogical(:,:,i);
%                 slicePrior = prior(:,:,i);
%                 [~, ~, ~, ~, posteriors] = fitGMMwithPrior(...
%                     sliceData(sliceMask), ...
%                     priorsFiltered, ...
%                     priorMap, ...
%                     means, ...
%                     variances, ...
%                     weights, ...
%                     1);
%                 
%                 % threshold at maximum probability
%                 segmented(sliceMask) = max(posteriors, [], 2);
%             end

        case  'mixturewithsupervoxel'
            positionsX = repmat(repmat(((1:origDataSize(1))' - 0.5), [1, origDataSize(2)]), [1, 1, origDataSize(3)]);
            positionsY = repmat(repmat(((1:origDataSize(2)) - 0.5), [origDataSize(1), 1]), [1, 1, origDataSize(3)]);
            positionsZ = reshape(ones(origDataSize(1)*origDataSize(2), 1) * ((1:origDataSize(3)) - 0.5), [origDataSize(1), origDataSize(2), origDataSize(3)]);
            
            origDataPositions = [positionsX(:), positionsY(:), positionsZ(:)];
            clear positionsX positionsY positionsZ 
            
            % one gmm iteration to obtain the posteriors for all voxels
            origDataPositions = origDataPositions .* repmat(origStepSize, [numTotalOrigVoxels, 1]);
            
            numelPerSlize = origDataSize(1)*origDataSize(2);
            for i = 1:origDataSize(3)

                tmpSliceData = log(correctedData(:,:,i));
                sliceData = [tmpSliceData(:), origDataPositions((1+(i-1)*numelPerSlize):(i*numelPerSlize),:)];
                sliceMask = origMaskLogical(:,:,i);


                [~, ~, ~, ~, posteriors] = fitGMMwithSV(...
                    sliceData(sliceMask, :), ...
                    means, ...
                    variances, ...
                    weights, ...
                    superVoxelTheta, ...
                    updateVariance, ...
                    updateSpatialVariance, ...
                    updateWeights, ...
                    sameVariance, ...
                    1);
                
                % threshold at maximum probability
                [~,idx] = max(posteriors, [], 2);
                tmpSlice = segmented(:,:,i);
                tmpSlice(sliceMask) = idx;
                segmented(:,:,i) = tmpSlice;
            end
            
        otherwise
            error('invalid Gaussian Mixture Model Fitting');
    end    
end

% end of iic
end


