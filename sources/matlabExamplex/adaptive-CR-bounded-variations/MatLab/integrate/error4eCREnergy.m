function val = error4eCREnergy(c4n,n4e,gradExact,uApprox)
    nrElems = size(n4e,1);
    s4e = computeS4e(n4e);
    gradUApprox = zeros(nrElems,2);
    
    % Compute gradient for x
    for elem = 1:size(n4e,1);
        grads = [c4n(n4e(elem,:),:)'; 1 1 1] \ [-2 0; 0 -2; 0 0];
        gradUApprox(elem,:) = uApprox(s4e(elem,:))' * grads([3 1 2],:);
    end
    
    val = sum(integrate(c4n,n4e, ...
        @(n4p, Gpts4p, Gpts4ref) (gradExact(Gpts4p) - gradUApprox).^2, 6),2);
end