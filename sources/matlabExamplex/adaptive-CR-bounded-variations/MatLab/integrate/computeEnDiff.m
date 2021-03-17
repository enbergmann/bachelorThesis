function enDiff = computeEnDiff(sigmaCR4e, W,...
                    DW, m, c4n, n4e, area4e, minEnergy, intF4e, WConjugate)
%% computeEnDiff - compute the difference between Ed* and E* of sigmadRT
%
% Input:     W              convex energy density
%            DW             derivative of W
%            c4n            coordinates for the nodes of the mesh
%            n4e            nodes for the elements of the mesh
%            minEnergy      minimum energy
%            intF4e         integral of F for each element
%            area4e         area of each element
%
% Optional   WConjugate     convex conjugate of W
%
% Output:    endiff         Ed*(sigmadRT) - E*(sigmadRT)

    nrElems = size(n4e,1);
    meanF4e = intF4e./repmat(area4e,1,m);
    mid4e = (c4n(n4e(:,1),:) + c4n(n4e(:,2),:) + c4n(n4e(:,3),:))/3;
    degWConjugate = 10; % degree for Gauss integration of the convex conjugate
    
    EnSigmadRT4e = zeros(nrElems,1);
    
    if nargin < 10
        WConjugate = @(x)computeWConjugate(x, W, DW, m); % if no convex conjugate is given
    end
    
    for j = 1:nrElems
        sigmaCR = sigmaCR4e(j,:);
        meanF = meanF4e(j,:);
        mid = mid4e(j,:);
        
        EnSigmadRT4e(j) = integrate(c4n,n4e(j,:),@(n4p,Gpts4p,Gpts4ref)...
            integrandWConjugate(Gpts4p,sigmaCR,WConjugate,m,meanF,mid),degWConjugate,area4e(j));
    end
    enDiff = minEnergy + sum(EnSigmadRT4e);
end

function val = computeWConjugate(x, W, DW, m)
%% computeWConjugate - compute the convex conjugate of W via quasi-Newton (BFGS)

    nrComps = size(x,1);
    val = zeros(nrComps,1);
    for j = 1:nrComps
        xStar = x(j,:);
        options = optimoptions('fminunc','Display','off',...
                    'Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
        problem.options = options;
        problem.x0 = zeros(1,2*m);
        problem.objective = @(y)objectiveFunctional(y, xStar, W, DW);
        problem.solver = 'fminunc';
        [~,min] = fminunc(problem);
        val(j) = -min;
    end
end

function [f,df] = objectiveFunctional(x, xStar, W, DW)
    f = - sum(x.*xStar,2) + W(x);
    df = - xStar + DW(x);
end

function val = sigmadRT(x, sigmaCR, m, meanF, mid)
%% sigmadRT - compute sigmadRT = sigmaCR - Pi0 f/2 (x - mid)

    val = zeros(size(x,1),2*m);
    for j = 1:size(x,1)
        for k = 1:m
            val(j,[2*k-1,2*k]) = sigmaCR([2*k-1,2*k]) - meanF(k)*(x - mid)/2;
        end
    end
end

function val = integrandWConjugate(x, sigmaCR, WConjugate, m, meanF, mid)
    y = sigmadRT(x, sigmaCR, m, meanF, mid);
    val = WConjugate(y);
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!