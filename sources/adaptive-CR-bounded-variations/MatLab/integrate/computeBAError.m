function [errF4e, mean4e] = computeBAError(c4n,n4e,f,degree,area4e,intF4e)
%% computeBAError - bestapproximation error in piecewise constants
% This functions computes the L2 error ||f-mean(f)||_L2(T)^2 of the
% bestapproximation of a given L2 function f in the space of the piecewise
% constant functions on the triangulation defined by c4n and n4e.
%
% COMMAND
%   errF4e=computeBAError(c4n,n4e,f,degree,area4e,intF4e)
%
% INPUT
%   c4n,n4e ... triangulation
%   f       ... L2 function
%   degree  ... polynomial degree for Gauss quadrature
%   area4e  ... area of the the triangles (OPTIONAL)
%   intF4e  ... piecewise constant approximation of f
%
% OUTPUT
%   errF4e  ... bestapproximation error
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

%% PROCEED INPUT
if nargin<5; area4e=computeArea4e(c4n,n4e); end
if nargin<6;
  % Compute the integral of f on each element
  intF4e=integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p)),degree);
end

%% INITIALIZATION
nrComps=length(f(c4n(n4e(1),:)));
oscDeg=2*degree;

%% COMPUTATION
% Compute integral mean
mean4e=intF4e./repmat(area4e,1,nrComps);
% bestapproximation error
errF4e=integrate(c4n,n4e,...
                 @(n4p,Gpts4p,Gpts4ref)sum((f(Gpts4p)-mean4e).^2,2),oscDeg);
end


