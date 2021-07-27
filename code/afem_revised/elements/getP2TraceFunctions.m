function [ baseQ ] = getP2TraceFunctions(  )
%GETP2TRACEFUNCTIONS Create P0 basis function.
%   GETP2TRACEFUNCTIONS() returns a structure array containing the
%   function handle of the four P0 basis function on the reference side.
%
% output:   baseQ       -   lagrange basis on the reference side for degree
%                           2 with uniform steps

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

baseQ = cell(3, 1);
baseQ{1} = @(x) 2*(x-0.5).*(x-1);
baseQ{2} = @(x) 4*x.*(1-x);
baseQ{3} = @(x) 2*x.*(x-0.5);
end

