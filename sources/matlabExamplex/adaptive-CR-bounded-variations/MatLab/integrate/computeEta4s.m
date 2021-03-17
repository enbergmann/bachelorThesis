function [eta4s,n4s] = computeEta4s(sigma4e,c4n,n4e,s4e,n4sDb,q)

    %% Initialisation
    n4s      = computeN4s(n4e);
    length4s = computeLength4s(c4n,n4s);
    
    %% Compute the L2-norm of the jumps and weigh them with length4s.
    % eta4sNormal = P0NormalJump(c4n,n4e,n4sDb,n4s,s4e,sigma4e);
    eta4sTangent = P0TangentJump(c4n,n4e,sigma4e);
    % eta4s = length4s.^2.*(eta4sNormal + eta4sTangent).^q;
    eta4s = length4s.^2.*eta4sTangent.^q;
end