function [exactErrorCR, exactErrordRT] ...
            = computeExactError(sigmaCR4e, DuExact, DW, c4n, n4e, r, t, intF4e, area4e)
%% computeExactError - compute exact error || sigma - sigmaCR ||_r/t and || sigma - sigmadRT ||_r/t
%
% Input:     sigmaCR4e      discrete solution sigmaCR on each element
%            DuExact        derivative of exact solution
%            DW             derivative of W
%            c4n            coordinates for nodes
%            n4e            nodes for elements
%            r, t           parameters
%            area4e         area of each element
%
% Output:    exactError     || sigma - sigmaCR ||_r/t
%            exactErrordRT  || sigma - sigmadRT ||_r/t

    m = size(intF4e,2);
    nrElems = size(n4e,1);
    meanF4e = intF4e./repmat(area4e,1,m);
    mid4e = (c4n(n4e(:,1),:) + c4n(n4e(:,2),:) + c4n(n4e(:,3),:))/3;
    
    % || sigma - sigmaCR ||_r/t
    exactError4e = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                    ((sum((DW(DuExact(Gpts4p))-sigmaCR4e).^2,2)).^(r/(2*t))),8,area4e);
    exactErrorCR = sum(exactError4e)^(t/r);
    
    % || sigma - sigmadRT ||_r/t
    exactErrordRT4e = zeros(1,nrElems);
    for j = 1:nrElems
        sigmaCR = sigmaCR4e(j,:);
        meanF = meanF4e(j,:);
        mid = mid4e(j,:);
        
        exactErrordRT4e(j) = integrate(c4n,n4e(j,:),@(n4p,Gpts4p,Gpts4ref)...
                               ((sum((DW(DuExact(Gpts4p))-sigmadRT(Gpts4p, ...
                                sigmaCR, m, meanF, mid)).^2,2)).^(r/(2*t))),8,area4e(j));
    end
    exactErrordRT = sum(exactErrordRT4e)^(t/r);
end

function val = sigmadRT(x, sigmaCR, m, meanF, mid)
%% sigmadRT - compute sigmadRT = sigmaCR - Pi0 f/2 x
    val = zeros(size(x,1),2*m);
    for j = 1:size(x,1)
        for k = 1:m
            val(j,[2*k-1,2*k]) = sigmaCR([2*k-1,2*k]) - meanF(k)*(x - mid)/2;
        end
    end
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!