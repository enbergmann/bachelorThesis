function [Phi, GradPhi,nCP, nIP , c4Base] = getP2ShapeFunctions()
%%GETP2SHAPEFUNCTION Creates P2 basis functions and their gradients.
%   GETP2SHAPEFUNCTION() returns a structure array containing the
%   function handle of the P2 basis functions on the reference triangle.
%
%   See also TRANSFORMAFFINE.
%
% output    Phi         - basis on the reference triangle for total degree 
%                         2; order: base functions in vertices,on first side,
%                         second side, third side, of inner points
%           GradPhi     - gradients of Phi (in the same order)
%           nCP         - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIP         - number of base functions in inner points

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

phi3 = @(x)( 1 - x(:,1) - x(:,2) ).* (2*( 1 - x(:,1) - x(:,2) )-1);
phi1 = @(x)( x(:,1) )             .* (2*x(:,1)-1);
phi2 = @(x)( x(:,2) )             .* (2*x(:,2)-1);

phi6 = @(x) 4*x(:,1)                 .*  x(:,2);
phi4 = @(x) 4*( 1 - x(:,1) - x(:,2) ).*  x(:,2);
phi5 = @(x) 4*( 1 - x(:,1) - x(:,2) ).*  x(:,1);


gradPhi3 = @(x) [-3 + 4*x(:,1)+4*x(:,2),-3 + 4*x(:,1)+4*x(:,2)];
gradPhi1 = @(x) [4*x(:,1)-1, zeros(size(x,1),1)];
gradPhi2 = @(x) [zeros(size(x,1),1), 4*x(:,2)-1];

gradPhi6 = @(x) 4* [x(:,2), x(:,1)];
gradPhi4 = @(x) 4* [-x(:,2),1-x(:,1)-2*x(:,2)];
gradPhi5 = @(x) 4* [1-2*x(:,1)-x(:,2), -x(:,1)];

Phi = {phi1, phi2, phi3, phi4, phi5, phi6};
GradPhi = {gradPhi1, gradPhi2, gradPhi3, gradPhi4, gradPhi5, gradPhi6};

nCP = 1; nIP = 0; c4Base= [0.5];
end