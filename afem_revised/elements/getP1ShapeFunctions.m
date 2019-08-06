function [Phi, GradPhi,nCP, nIP , c4Base] = getP1ShapeFunctions()
%%GETP1SHAPEFUNCTIONS Create P1 basis functions and their gradients.
%   GETP1SHAPEFUNCTIONS() returns a structure array containing the
%   function handle of the three P1 basis function on the reference triangle.
%
%   [Phi, GradPhi] = GETP1SHAPEFUNCTIONS() also returns the gradients of
%   the basis functions.
%
%   See also TRANSFORMAFFINE.
%
% output    Phi         - basis on the reference triangle for total degree 
%                         1; order: base functions in vertices,on first side,
%                         second side, third side, of inner points
%           GradPhi     - gradients of Phi (in the same order)
%           nCP         - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIP         - number of base functions in inner points

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

phi3 = @(x)( 1 - x(:,1) - x(:,2) );
phi1 = @(x)( x(:,1) );
phi2 = @(x)( x(:,2) );

gradPhi3 = @(x)repmat([-1, -1], size(x, 1), 1);
gradPhi1 = @(x)repmat([1, 0], size(x, 1), 1);
gradPhi2 = @(x)repmat([0, 1], size(x, 1), 1);

Phi = {phi1, phi2, phi3};
GradPhi = {gradPhi1, gradPhi2, gradPhi3};

nCP = 1; nIP = 0; c4Base = [];
end

