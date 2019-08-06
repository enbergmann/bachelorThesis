function [ baseQ ] = getP0TraceFunctions()
%GETP0TRACEFUNCTIONS Create P0 basis function.
%   GETP0TRACEFUNCTIONS() returns a structure array containing the
%   function handle of the P0 basis function on the reference side.
%
% output:   baseQ       -   lagrange basis on the reference side for degree
%                           0 with uniform steps

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

baseQ = cell(1, 1);
baseQ{1} = @(x) ones(size(x,1),1);

end

