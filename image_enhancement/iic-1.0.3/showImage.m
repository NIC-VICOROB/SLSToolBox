function  showImage(data, range, location)
%
%
% Original Author: Koen Van Leemput (kvle)
% Modified by Christian Thode Larsen, christian.thode.larsen@gmail.com

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

if (nargin < 2)
  range = double([min(data(:)), max(data(:))]);
end

[Nx, Ny, Nz] = size(data);

if (nargin < 3)
    location = [round(Nx * 0.5), round(Ny * 0.5), round(Nz * 0.5)];
end

xySlice = squeeze(data(:, :, location(3), :));
xzSlice = squeeze(data(:, location(2), :, :));
yzSlice = squeeze(data(location(1), :, :, :));


if (ndims(data) < 4)
  patchedSlices = [xySlice, xzSlice; yzSlice', zeros(Nz) + range(1)];
  patchedSlices = (single(patchedSlices) - range(1)) / (range(2) - range(1));
  imshow(patchedSlices',  'InitialMagnification', 'fit');
else
  % RGB data
  patchedSlices = [xySlice, xzSlice; permute(yzSlice, [2, 1, 3]), zeros(Nz, Nz, 3)];  
  imshow(permute(patchedSlices, [2, 1, 3]), 'InitialMagnification', 'fit');
end
