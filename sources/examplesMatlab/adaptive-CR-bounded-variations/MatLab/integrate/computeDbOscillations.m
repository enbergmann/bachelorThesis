function [osc4sDb,mean4sDb,oscDb4e]=...
            computeDbOscillations(c4n,n4e,n4sDb,Dg,degree)
%% computeDbOscillations - oscillations of arc-length derivative on the boundary
% This function returns oscillations of the arc-length derivative of some
% given data on the Dirichlet boundary.
%
% COMMAND
%   [osc4sDb,mean4sDb,oscDb4e]=computeDbOscillations(c4n,n4e,n4sDb,Dg,degree)
%
% INPUT
%   c4n,n4e   ... triangulation
%   n4sDb     ... nodes for the Dirichlet boundary sides
%   Dg        ... gradient of the boundary data
%                 (gradients of all components in one row)
%   degree    ... polynomial degree for Gauss quadrature
%
% OUTPUT
%   osc4sDb   ... oscillations on the boundary
%   mean4sDb  ... integral mean of the arc-length derivative
%   osc4e     ... elementwise accumulation of the oscillations
%
% See also integrate.

% Copyright 2014 Philipp Bringmann
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
nrComps=length(Dg(c4n(n4sDb(1),:)))/2;
degreeOsc=2*degree;
tangent4sDb=computeTangent4s(c4n,n4sDb);
length4sDb=computeLength4s(c4n,n4sDb);

%% COMPUTATION
% integral mean
mean4sDb=...
  integrate(c4n,n4sDb,@(n4p,Gpts4p,Gpts4ref)Dg(Gpts4p).*...
            repmat(tangent4sDb,1,nrComps),degree)./...
  repmat(length4sDb,1,2*nrComps);
% oscillations
osc4sDb=...
  integrate(c4n,n4sDb,@(n4p,Gpts4p,Gpts4ref)...
            sum(((Dg(Gpts4p).*repmat(tangent4sDb,1,nrComps)-mean4sDb).^2),2),...
            degreeOsc).*length4sDb;
% accumulate integral mean
mean4sDb=mean4sDb(:,1:2:end)+mean4sDb(:,2:2:end);

%% ELEMENTWISE ACCUMULATION
if nargout>2
  s4n=computeS4n(n4e);
  e4s=computeE4s(n4e);
  DbSides=rowaddr(s4n,n4sDb(:,1),n4sDb(:,2));
  e4sDb=e4s(DbSides,1);
  oscDb4e=accumarray(e4sDb,osc4sDb,[size(n4e,1),1]);
end

end
