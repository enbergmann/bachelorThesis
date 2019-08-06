function [GradPhi4sPlus, GradPhi4sMinus] = ...
  transformJumpsAffine(c4n, n4e, GradPhi)
%% TRANSFORMJUMPSAFFINE Affine transformation of gradient jumps on edges.
%   [GradPhi4sPlus, GradPhi4sMinus] = TRANSFORMJUMPSAFFINE(c4n, n4e, GradPhi)
%   returns the affine transformation of the gradient function handles in
%   the cell array GradPhi on the two adjacent triangles of each edge. The
%   arrays c4n and n4e describe the triangulation. The function handles in
%   the cell arrays GradPhi4sPlus and GradPhi4sMinus expect 1D Gauss points
%   as arguments. Each function handle in GradPhi4sPlus returns a
%   (no of sides) x 2 array and each function handle in GradPhi4sMinus a
%   (no of interior sides) x 2 array.
%
%   See also TRANSFORMAFFINE.

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
nFunctions = length(GradPhi);
n4s = computeN4s(n4e);
e4s = computeE4s(n4e);
nSides = size(n4s, 1);
c4e = [c4n(n4e(:,1),:), c4n(n4e(:,2),:), c4n(n4e(:,3),:)];
% jump data
ePlus = e4s(:,1);
isInterior = e4s(:,2)~=0;
eMinus = e4s(isInterior,2);

%% CONSTANTS
% data structure for reference triangle
c4sRef = [1 0 0 1; 0 1 0 0; 0 0 1 0];
n4sRef = [2 3; 3 1; 1 2];

%% COMPUTE EDGE COORDINATES ON REFERENCE TRIANGLE
% computes coordinates on the reference triangle of the two nodes of the
% given edge
c4sP = zeros(nSides, 4);  % PLUS
c4sM = zeros(nnz(isInterior), 4);  % MINUS
for j = 1:3
  IndPlus = all(n4e(ePlus,n4sRef(j,:))==n4s, 2) ...
            | all(n4e(ePlus,n4sRef(j,:))==n4s(:,[2 1]), 2);
  c4sP(IndPlus,:) = repmat(c4sRef(j,:), nnz(IndPlus), 1);
  IndMinus = ...
    all(n4e(eMinus,n4sRef(j,:))==n4s(isInterior,:), 2) ...
    | all(n4e(eMinus,n4sRef(j,:))==n4s(isInterior,[2 1]), 2);
  c4sM(IndMinus,:) = repmat(c4sRef(j,:), nnz(IndMinus), 1);
end

%% AFFINE TRANSFORMATION AND INVERSION
TRAFO4e = [c4e(:,[3 5])-c4e(:,[1 1]), c4e(:,[4 6])-c4e(:,[2 2])];
DET4e = TRAFO4e(:,1).*TRAFO4e(:,4) - TRAFO4e(:,2).*TRAFO4e(:,3);
INVTRAFO4ePlus = [ TRAFO4e(ePlus,4), -TRAFO4e(ePlus,2),...
                  -TRAFO4e(ePlus,3),  TRAFO4e(ePlus,1)]...
                 ./ repmat(DET4e(ePlus), 1, 4);
INVTRAFO4eMinus = [ TRAFO4e(eMinus,4), -TRAFO4e(eMinus,2),...
                   -TRAFO4e(eMinus,3),  TRAFO4e(eMinus,1)]...
                  ./ repmat(DET4e(eMinus), 1, 4);

%% TRANSFORM GRADIENTS
% transform 1D Gauss point to 2D Gauss point on reference triangle and
% transform gradients to the specific triangulation
GradPhi4sPlus = cell(nFunctions, 1);
GradPhi4sMinus = cell(nFunctions, 1);
for j = 1:nFunctions
  GradPhi4sPlus{j} = @(x) ...
    [sum(GradPhi{j}(c4sP(:,[1 2])+x*(c4sP(:,[3 4])-c4sP(:,[1 2]))) ...
         .*INVTRAFO4ePlus(:,[1 3]), 2), ...
     sum(GradPhi{j}(c4sP(:,[1 2])+x*(c4sP(:,[3 4])-c4sP(:,[1 2]))) ...
         .*INVTRAFO4ePlus(:,[2 4]), 2)];
  GradPhi4sMinus{j} = @(x) ...
    [sum(GradPhi{j}(c4sM(:,[1 2])+x*(c4sM(:,[3 4])-c4sM(:,[1 2]))) ...
         .*INVTRAFO4eMinus(:,[1 3]), 2), ...
     sum(GradPhi{j}(c4sM(:,[1 2])+x*(c4sM(:,[3 4])-c4sM(:,[1 2]))) ...
         .*INVTRAFO4eMinus(:,[2 4]), 2)];
end


end
