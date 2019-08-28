function eta4e = estimateError4e(u,f,c4n,n4e,n4sDb,n4sNb,alpha,delta)
% TODO name should contain CR 
  area4e = computeArea4e(c4n,n4e);
  n4s = computeN4s(n4e);
  s4e = computeS4e(n4e);
  nrElems = size(n4e,1);
  mid4e = computeMid4e(c4n,n4e);
  length4s = computeLength4s(c4n,n4s);
  nodeValues4e = computeNodeValues4e(s4e,u);
  s4n = computeS4n(n4e,n4s);
  e4s = computeE4s(n4e);
    % nodeValues4e(j) = (u(P1),u(P2),u(P3)) wrt. T_j
  absNodeJumps4s = computeAbsNodeJumps4s(n4e,e4s,nodeValues4e);

  n=2; %?
  [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,20,area4e);
    % temp_j(elem) = \int_T psi_j *f dx
  termF2 = integrate(@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p).^2),c4n,n4e,20,area4e);
    % temp(elem) = ||f||^2_{L^2(T)}

  termMixed = u(s4e(:,1)).*temp1+u(s4e(:,2)).*temp2+u(s4e(:,3)).*temp3;
  termU = 1/3*area4e.*sum(u(s4e).^2,2);
  termJumps = 1/4*(...
              length4s(s4e(:,1)).*sum(absNodeJumps4s(s4e(:,1),:),2)...
              + length4s(s4e(:,2)).*sum(absNodeJumps4s(s4e(:,2),:),2)...
              + length4s(s4e(:,3)).*sum(absNodeJumps4s(s4e(:,3),:),2)); 


  %eta4e = area4e.^(2/n).*(termF-2*alpha*(...
  %  u(s4e(:,1))'*temp1+u(s4e(:,2))'*temp2+u(s4e(:,3))'*temp3)...
  %  + alpha^2/3*area4e.*sum(u(s4e),2))... % Volumenanteil
  %  + 1/4*area4e.^(delta/n).*(...
  %   length4s(s4e(:,1))'*sum(absNodeJumps4s(s4e(:,1)),2)...
  %  +length4s(s4e(:,2))'*sum(absNodeJumps4s(s4e(:,2)),2)...
  %  +length4s(s4e(:,3))'*sum(absNodeJumps4s(s4e(:,3)),2)); %Sprunganteil

  eta4e = area4e.^(2/n).*(termF2-2*alpha*termMixed+alpha^2*termU)...
          + area4e.^(delta/n).*termJumps;
end
