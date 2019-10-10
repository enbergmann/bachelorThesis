function eta4e = estimateErrorCR4e(params, currData, u)
 %TODO some of those should prob. be in currData
  alpha4Estimate = params.alpha4Estimate;
  beta4Estimate = params.beta4Estimate;
  f = params.f;

  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  n4s = currData.n4s;
  s4e = currData.s4e;
  nrElems = currData.nrElems;
  mid4e = computeMid4e(c4n,n4e);
  length4s = currData.length4s;
  s4n = computeS4n(n4e,n4s);
  e4s = computeE4s(n4e);
    % nodeValues4e(j) = (u(P1),u(P2),u(P3)) wrt. T_j

    % TODO prob rewrite some of this stuff, at least documentation
  nodeValues4e = computeNodeValues4e(s4e, u);
  absNodeJumps4s = computeAbsNodeJumps4s(n4e, e4s, nodeValues4e);

  n=2; %?
  int1RHS4e = currData.int1RHS4e;
  int2RHS4e = currData.int2RHS4e;
  int3RHS4e = currData.int3RHS4e;
    % temp_j(elem) = \int_T psi_j *f dx
  termF2 = integrate(@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p).^2),c4n,n4e,20,area4e);
    % temp(elem) = ||f||^2_{L^2(T)}

  termMixed = u(s4e(:,1)).*int1RHS4e + u(s4e(:,2)).*int2RHS4e ... 
    + u(s4e(:,3)).*int3RHS4e;
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

  eta4e = area4e.^(2/n).*(termF2 - 2*alpha4Estimate*termMixed + ...
    alpha4Estimate^2*termU)...
          + area4e.^(beta4Estimate/n).*termJumps;
end
