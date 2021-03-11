function error4e = error4eP1Energy(c4n,n4e,gradExact,uApprox)
%% error4eP1Energy - compute the energy error on elements
%  Input:  c4n,n4e     - mesh
%          gradExact   - the exact gradient
%          uApprox     - AFEM P1 solution
%
%  Output: error4e     - the exact squared energy error on each element


    %% Compute grad U
    nrElems = size(n4e,1);
    gradU = zeros(nrElems,2);

    for elem = 1 : nrElems
        grads = [1,1,1;c4n(n4e(elem,:),:)']\[0,0;eye(2)];
        gradU(elem,:) = uApprox(n4e(elem,:))' * grads;
    end

    %% Compute the error
    
    error4e = integrate(c4n,n4e,@(n4p, Gpts4p, Gpts4ref) (...
	        sum((gradExact(Gpts4p) - gradU ).^2,2)),6);                    
    
    
end        

