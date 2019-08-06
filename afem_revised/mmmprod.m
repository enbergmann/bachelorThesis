function val = mmmprod(A, B, C)
%%MMMPROD Computes triple matrix products
%   MMMPROD(A, B, C) returns the triple product of MxM matrices
%   A' * B * C where A, B, and C are two-dimensional N x (M^2) arrays, e.g.
%   A = [A11, A12, A21, A22]. The output dimension is determined by B.
%
%   See also VMVPROD.

% Copyright 2016 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

oldDim = size(B);

if ismatrix(A)
  A = reshape(A, [size(A, 1), sqrt(size(A, 2)), sqrt(size(A, 2))]);
  A = permute(A, [1 3 2]);
end
if ismatrix(B)
  B = reshape(B, [size(B, 1), sqrt(size(B, 2)), sqrt(size(B, 2))]);
  B = permute(B, [1 3 2]);
end
if ismatrix(C)
  C = reshape(C, [size(C, 1), sqrt(size(C, 2)), sqrt(size(C, 2))]);
  C = permute(C, [1 3 2]);
end

val = bsxfun(@times, permute(A,[1 2 4 3]), permute(B,[1 4 3 2]));
val = sum(val, 4);
val = bsxfun(@times, permute(val,[1 2 4 3]), permute(C,[1 4 3 2]));
val = sum(val, 4);

val = reshape(val, oldDim);

end
