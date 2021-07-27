function val = vmvprod(A,u,v)
%%VMVPROD Computes triple matrix products
%   VMVPROD(A, u, v) returns the vector-matrix-vector linewise product
%   u' * A * v, where A is an M x M x N tensor and u and v are N x M
%   arrays.
%
%   See also MMMPROD.

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

val = bsxfun(@times, permute(A, [3 1 2]), v);
val = squeeze(sum(val, 2));
val = sum(val.*u, 2);

end
