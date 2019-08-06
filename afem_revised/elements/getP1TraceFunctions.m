function [ baseQ ] = getP1TraceFunctions()
%GETP1TRACEFUNCTIONS Create P1 basis function.
%   GETP1TRACEFUNCTIONS() returns a structure array containing the
%   function handle of the two P1 basis function on the reference side.
%
% output:   baseQ       -   lagrange basis on the reference side for degree
%                           1 with uniform steps

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a
                 
baseQ = cell(2, 1);
baseQ{1} = @(x) 1-x;
baseQ{2} = @(x) x;

end

