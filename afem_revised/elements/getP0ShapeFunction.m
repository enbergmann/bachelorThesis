function [Phi, GradPhi, nCP, nIP , c4Base] = getP0ShapeFunction()
%%GETP0SHAPEFUNCTION Create P0 basis function.
%   Phi = GETP0SHAPEFUNCTION() returns a structure array containing the
%   function handle of the one P0 basis function on the reference triangle.
%
%   See also TRANSFORMAFFINE.
%
% output    Phi         - basis on the reference triangle for total degree 
%                         0; order: base functions in vertices,on first side,
%                         second side, third side, of inner points
%           GradPhi     - gradients of Phi (in the same order)
%           nCP         - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIP         - number of base functions in inner points

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

phi1 = @(x) ones(size(x, 1), 1);
gradPhi1 = @(x) zeros(size(x,1),1);

Phi = {phi1}; GradPhi = {gradPhi1};

nCP =0; nIP = 1; c4Base = [];
end

