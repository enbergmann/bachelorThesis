function [uCR, Energy, step] = solvePrimalDualFormulation(uCR, Lambda4e, alpha, epsilon, tau, ...
                                STIMACR, MASSCR, intFCR14e, intFCR24e, intFCR34e,...
                                gradsCR4e, s4e, dof, nrSides, area4e)
    
    nrElems = size(s4e,1);                        
    A = STIMACR/tau + alpha*MASSCR;
    
    vCR = zeros(nrSides,1);
    
    DuCR4e = computeGradientCR(uCR, gradsCR4e, s4e);
    step = 1;
    
    while true
        DvCR4e = computeGradientCR(vCR, gradsCR4e, s4e);
        
        % utilde = u_{k-1} + tau*v_{k-1}
        DuTilde = DuCR4e + tau*DvCR4e;
        tmp = Lambda4e + tau*DuTilde;
        % Lambda_{k} = (Lambda_{k-1} + tau*Grad utilde)/max{1, |Lambda_{k-1} + tau*Grad utilde|}
        Lambda4e = tmp./repmat(max(1,sqrt(sum(tmp.^2, 2))), 1, 2);
        
        RHS = zeros(nrSides, 1);
        for elem = 1 : nrElems
            sides = s4e(elem,:);
            area = area4e(elem);
            
            % RHS = (Grad u_{k-1}, Grad phi_CR)/tau + (f, phi_CR) - (Lambda_k, Grad phi_CR)
            RHS(sides) = RHS(sides)...
                + area.*sum(repmat(DuCR4e(elem,:),3,1).*gradsCR4e(:,:,elem), 2)/tau...
                + [intFCR14e(elem); intFCR24e(elem); intFCR34e(elem)] ...
                - area.*sum(repmat(Lambda4e(elem,:),3,1).*gradsCR4e(:,:,elem), 2);
            
        end
        
        uCRnew = zeros(nrSides, 1);
        uCRnew(dof) = A(dof,dof)\RHS(dof);

        vCR = (uCRnew - uCR)/tau;        

        % stop criteria: ||Grad(u_k - u_{k-1})/tau|| < epsilon
        % ||Grad(u_k - u_{k-1})/tau||^2 = sum over T of |T|*|Grad(u_k - u_{k-1})|^2/tau^2
        DuCRnew4e = computeGradientCR(uCRnew, gradsCR4e, s4e);

        tmp = (DuCRnew4e - DuCR4e)/tau;
        stop = sqrt(sum(area4e.*sum(tmp.^2, 2)));
        
        uCR = uCRnew;
        DuCR4e = DuCRnew4e;
        
        if stop < epsilon
            break;
        end
        
        step = step + 1;
    end
    
    % E_nc(uCR) = alpha/2*||uCR||^2 + int_Omega |D uCR| dx - int_Omega f uCR dx
    Energy4e = area4e.*sqrt(sum(DuCR4e.^2,2)) + alpha/6*area4e.*sum(uCR(s4e).^2, 2)...
                    - sum(uCR(s4e).*[intFCR14e, intFCR24e, intFCR34e],2);
    Energy = sum(Energy4e);
end