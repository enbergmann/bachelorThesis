function [osc4e ,mean4e] = oscillations(c4n,n4e,f,degree)
%% oscillcations - oscillations of a function on a mesh
%  Input:    c4n,n4e - mesh
%            f       - R^2->R; input: points; output: values
%            degree  - accuracy of integration
%  Output:   osc4e   - vector of oscillations squared for each element
%	     mean4e  - vector of integral means of f for each element
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

%% Initialisation
    meanDeg = degree;
    oscDeg  = 2*degree;
    area4e  = computeArea4e(c4n,n4e);
%% Compute the integral mean of f on each element.
    mean4e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)(f(Gpts4p)),meanDeg)...
             ./area4e;
%% Compute oscillations: locally on each T - osc4e = ||f-mean(f)||_L2(T) / |T|
    osc4e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)((f(Gpts4p)-mean4e).^2),...
    	     oscDeg);
    osc4e = osc4e .* area4e;
end


