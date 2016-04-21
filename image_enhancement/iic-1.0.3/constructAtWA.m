function AtWA = constructAtWA(weights, B2, indexMapB2)
% Description: Computes the Hessian using separability to speed up things
%
% Call:     hessian = computeHessian(weights, B2, B2Index)
%
% Input:    weights:    A 3D volume that contains voxel weights. If weights are equal, pass mask matrix.   
%           B2:         A cell containing three 2D matrices of size NxM^2 for a given direction, where N is the number of
%                       voxels and M^2 is the squared number of basis functions
%           indexMapB2: A precomputed index map which is necessary due to our separability speed-up.
%                       Obtained from computeSeparableBasisFunctions.
%
% Output:   AtWA:    AtWA (the hessian) which has size P*Q*R x P*Q*R, where P, Q and R is the number of basis functions in each
%                       direction

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

[Nx, Nbx] = size(B2{1});
[Ny, Nby] = size(B2{2});
[Nz, Nbz] = size(B2{3});

NbxyzSQRT = sqrt(Nbx * Nby * Nbz);

% then do AWA super fast and with sick indicing
% x direction
tmp = permute(reshape(B2{1}' * reshape(weights, [Nx, Ny * Nz]), [Nbx, Ny, Nz]), [2, 1, 3]);

% y direction
tmp = permute(reshape(B2{2}' * reshape(tmp, [Ny, Nbx * Nz]), [Nby, Nbx, Nz]), [2, 1, 3]);

% z direction
tmp = reshape(reshape(tmp, [Nbx * Nby, Nz]) * B2{3}, [Nbx, Nby, Nbz]);

AtWA = zeros(NbxyzSQRT, NbxyzSQRT);
AtWA(indexMapB2) = tmp(:);

end


