function [Phi4e, GradPhi4e] = transformAffine(c4n, n4e, Phi, GradPhi)
%%TRANSFORMAFFINE Affine transformation of functions on triangles.
%   [Phi4e] = TRANSFORMAFFINE(c4n, n4e, Phi) returns the affine
%   transformation of all functions in the cell array Phi to each
%   triangle of the triangulation given by c4n and n4e. The function
%   handles in the cell array Phi4e expects arguments on the reference
%   triangle and yield the corresponding function value on each triangle.
%
%   [Phi4e, GradPhi4e] = TRANSFORMAFFINE(c4n, n4e, Phi, GradPhi) also
%   transforms the gradients of the functions.
%
%   See also TRANSFORMJUMPSAFFINE.

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

%% INITIALIZATION
nFunctions = length(Phi);
nElem = size(n4e, 1);
c4e = [c4n(n4e(:,3),:), c4n(n4e(:,1),:), c4n(n4e(:,2),:)]; % or 1 2 3? 3 1 2

%% AFFINE TRANSFORMATION AND INVERSION
TRAFO4e = [c4e(:,[3 5])-c4e(:,[1 1]), c4e(:,[4 6])-c4e(:,[2 2])];
DET4e = TRAFO4e(:,1).*TRAFO4e(:,4) - TRAFO4e(:,2).*TRAFO4e(:,3);
INVTRAFO4e = [TRAFO4e(:,4), -TRAFO4e(:,2), -TRAFO4e(:,3), TRAFO4e(:,1)]...
             ./ repmat(DET4e, 1, 4);

%% TRANSFORM SHAPE FUNCTIONS
Phi4e = cell(nFunctions, 1);
%temp = INVTRAFO4e(:,[1 2]) + INVTRAFO4e(:,[3 4]);
for j = 1:nFunctions
  Phi4e{j} = @(x) repmat(Phi{j}(x), nElem, 1);
  %Phi4e{j} = @(x) Phi{j}(kron(x(:,1),temp(:,1))+kron(x(:,2),temp(:,2)));
end

%% TRANSFORM GRADIENTS
if nargin > 3
  GradPhi4e = cell(nFunctions, 1);
  for j = 1:nFunctions
    GradPhi4e{j} = ...
      @(x)[sum(repmat(GradPhi{j}(x), nElem, 1).*INVTRAFO4e(:,[1 3]), 2),...
           sum(repmat(GradPhi{j}(x), nElem, 1).*INVTRAFO4e(:,[2 4]), 2)];
  end
elseif nargout > 1
  error('Missing GradPhi input argument');
end

end
