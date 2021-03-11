function DuCR4e = computeGradientCR(uCR, gradsCR4e, s4e)
    nrElems = size(s4e,1);
    DuCR4e = zeros(nrElems,2);
    
    for elem = 1:nrElems
        sides = s4e(elem,:);
        DuCR4e(elem,:) = uCR(sides(1))*gradsCR4e(1, :, elem) + ...
            uCR(sides(2))*gradsCR4e(2, :, elem) + uCR(sides(3))*gradsCR4e(3, :, elem);
    end
end