function [eta, eta4e] =...
  estimatePoissonP1Residual(f, u4e, c4n, n4e, n4sDb, degreeF)

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
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
area4e = computeArea4e(c4n, n4e);
length4s = computeLength4s(c4n, n4s);
normal4s = computeNormal4s(c4n, n4s);

%% GET LOCAL SHAPE FUNCTIONS AND GRADIENTS
% P1
[~, GradPhi] = getP1ShapeFunctions();

%% TRANSFORM GRADIENT JUMPS ON EACH EDGE
[GradPhi4sPlus, GradPhi4sMinus] = ...
  transformJumpsAffine(c4n, n4e, GradPhi);

%% VOLUME CONTRIBUTION
INT = @(n4p, Gpts4p, Gpts4ref) f(Gpts4p).^2;
contrVol4e = area4e .* integrate(INT, c4n, n4e, 2*degreeF, area4e);

%% EDGE CONTRIBUTION
INT = @(n4p, Gpts4p, Gpts4ref) ...
      normalJump(Gpts4ref, GradPhi4sPlus, GradPhi4sMinus,...
                 u4e, e4s, normal4s).^2;
contrEdge4s = integrate(INT, c4n, n4s, 2, length4s);
contrEdge4e = sqrt(area4e) .* sum(contrEdge4s(s4e), 2);

%% ADD UP ALL CONTRIBUTIONS
eta4e = contrVol4e + contrEdge4e;
eta = sqrt(sum(eta4e));

end



function normalJump4s =...
  normalJump(Gpts4ref, GradPhi4sPlus, GradPhi4sMinus, u4e, e4s, normal4s)

%% INITIALIZATION
% jump data
ePlus = e4s(:,1);
isInterior = e4s(:,2)~=0;
eMinus = e4s(isInterior,2);

%% EVALUATE GRADIENTS AT GAUSS POINTS
GradPlus = [GradPhi4sPlus{1}(Gpts4ref), ...
            GradPhi4sPlus{2}(Gpts4ref), ...
            GradPhi4sPlus{3}(Gpts4ref)];
GradMinus = [GradPhi4sMinus{1}(Gpts4ref), ...
             GradPhi4sMinus{2}(Gpts4ref), ...
             GradPhi4sMinus{3}(Gpts4ref)];

%% COMPUTE NORMAL JUMPS
jump4s = [sum(u4e(ePlus,:).*GradPlus(:,[1 3 5]), 2),...
          sum(u4e(ePlus,:).*GradPlus(:,[2 4 6]), 2)];
jump4s(isInterior,:) = ...
  jump4s(isInterior,:) ...
  - [sum(u4e(eMinus,:).*GradMinus(:,[1 3 5]), 2),...
     sum(u4e(eMinus,:).*GradMinus(:,[2 4 6]), 2)];
jump4s(~isInterior,:) = 0;  % no contribution on boundary edges
normalJump4s =  sum(jump4s.*normal4s, 2);

end
