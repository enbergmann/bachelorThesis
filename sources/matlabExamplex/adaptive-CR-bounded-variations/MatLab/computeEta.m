function [etaSq, etaSq4e, volError4e, jumpError4e] ...
    = computeEta(uCR, alpha, beta, area4e, length4s, n4e, e4s, n4s,...
                 s4e, intFSq4e, intFCR14e, intFCR24e, intFCR34e)
    
    nrElems = size(area4e, 1);
    etaSq4e = zeros(nrElems, 1);
    volError4e = zeros(nrElems, 1);
    jumpError4e = zeros(nrElems, 1);
    nrSides = size(e4s, 1);
    
    uCR4e = zeros(nrElems, 3);
    uCRJump4s = zeros(nrSides, 2);
    
    for elem = 1:nrElems
        sides = s4e(elem, :);
        uloc = uCR(sides);
        uCR4e(elem, :) = uloc([1 2 3]) + uloc([3 1 2]) - uloc([2 3 1]);
    end
    
    for side = 1:nrSides
        Tp = e4s(side, 1);
        Tm = e4s(side, 2);
        node1 = n4s(side, 1);
        node2 = n4s(side, 2);
        
        ind1 = (n4e(Tp,:) == node1);
        ind2 = (n4e(Tp,:) == node2);
        
        uCRTp1 = uCR4e(Tp, ind1);
        uCRTp2 = uCR4e(Tp, ind2);
        
        if Tm == 0
            uCRTm1 = 0;
            uCRTm2 = 0;
        else
            ind1 = (n4e(Tm,:) == node1);
            ind2 = (n4e(Tm,:) == node2);
            uCRTm1 = uCR4e(Tm, ind1);
            uCRTm2 = uCR4e(Tm, ind2);
        end
        
        uCRJump4s(side,:) = [uCRTp1 - uCRTm1, uCRTp2 - uCRTm2];
    end
    
    for elem = 1:nrElems
        sides = s4e(elem,:);
        uCRloc = uCR(sides);
        area = area4e(elem);
        
        % first term |T|*||f - alpha*u||^2
        volError4e(elem) = area*intFSq4e(elem) ...
            + alpha^2*area^2*sum(uCRloc.^2)/3 ...
            - 2*alpha*area*sum(uCRloc.*[intFCR14e(elem); intFCR24e(elem); intFCR34e(elem)]);
        % second term |T|^(beta/2)*sum_{F in Fcal(T)} ||[uCR]_F||_L1(F)
        % ||[uCR]_F||_L1(F) = |F|(|[uCR(P1)]_F| + |[uCR(P2)]_F|)/4 with F = conv{P1, P2}
        jumpError4e(elem) = area^(beta/2)*...
                            sum(length4s(sides).*sum(abs(uCRJump4s(sides,:)),2))/4;
        
        etaSq4e(elem) = volError4e(elem) + jumpError4e(elem);
    end
    
    etaSq = sum(etaSq4e);
end