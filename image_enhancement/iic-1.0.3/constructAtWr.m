function AtWr = constructAtWr(data, weights, B)
% Author: Christian Thode Larsen (cthla), christian.thode.larsen@gmail.com

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

[Nx, NumBFuncsX] = size(B{1});
[Ny, NumBFuncsY] = size(B{2});
[Nz, NumBFuncsZ] = size(B{3});

% Calculate coefficients using separability of basis functions to speed things up

% x direction
tmp = permute(reshape(B{1}' * reshape(data .* weights, [Nx, Ny * Nz]), [NumBFuncsX, Ny, Nz]), [2, 1, 3]);

% y direction
tmp = permute(reshape(B{2}' * reshape(tmp, [Ny, NumBFuncsX*Nz]), [NumBFuncsY, NumBFuncsX, Nz]), [2, 1, 3]);

% z direction
AtWr = reshape(reshape(tmp, [NumBFuncsX * NumBFuncsY, Nz]) * B{3}, [NumBFuncsX, NumBFuncsY, NumBFuncsZ]);

AtWr = AtWr(:);
 % cthla: I'm abusing the fact that the smoothing operation is really the same, but for the weights, the hessian needs to be recomputed.
 % No need to do that here though, so hessian is provided as R   
% if(isempty(Rtranspose))
%     C(:) = R \ C(:);
% else
%     C(:) = R \ (Rtranspose \ C(:));
% end



