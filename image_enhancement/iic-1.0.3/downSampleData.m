function [downSampledData, start, stepSize, idx] = downSampleData(data, origStepSize, desiredStepSize, downSamplingType)

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

if nargin < 4
    downSamplingType = 'noInterp';
end

downSamplingType = lower(downSamplingType);

d0 = size(data);

switch downSamplingType
    case 'nointerp'
        % (From SPM): Discrete number of voxels to skip when downsampling
        sk = max([1, 1, 1], round(desiredStepSize * [1, 1, 1] ./ origStepSize));
        
        % downsample
        downSampledData = data(1:sk(1):d0(1), 1:sk(2):d0(2), 1:sk(3):d0(3));
        
        % voxels that are kept in original image
        idx = zeros(size(data));
        idx(1:sk(1):d0(1), 1:sk(2):d0(2), 1:sk(3):d0(3)) = 1;
        idx = idx > 0;
        
        % First voxel starts in 0,0,0 (center)
        start = [0,0,0];
        
        % Even if we specify a desiredStepSize, we're not guaranteed to get it if we do not interpolate
        stepSize = sk .* origStepSize;
        
    case 'nearestneighbor'
        d = ceil((d0-1)/desiredStepSize)+1;
        
        % The resampling assumes that we start and end at the outermost voxels, effectively moving the first and last
        % center to center -/+ 0.5*step respectively
        % round gives us the index that's closest to the voxel we're sampling for
        idxX = round(linspace(1,d0(1),d(1)));
        idxY = round(linspace(1,d0(2),d(2)));
        idxZ = round(linspace(1,d0(3),d(3)));              
        downSampledData = data(idxX, idxY, idxZ);
        
        % First voxel starts in 0,0,0 (center)
        start = [0,0,0];
        
        % This is actually not that precise...
        stepSize = desiredStepSize * [1,1,1];
    otherwise
        error('downsampling type not defined.');
end

end