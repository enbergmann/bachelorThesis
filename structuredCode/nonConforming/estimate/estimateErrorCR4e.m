function [eta4e, mu4e, xi4e] = estimateErrorCR4e(params, currData, output)
  %TODO interface documentation

  beta4Estimate = params.beta4Estimate;

  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  length4s = currData.length4s;
  e4s = currData.e4s;

  u = output.u;
  normOfDifference4e = output.normOfDifference4e;

  n=2; % TODO should probably in params, think about a suiting name and do it
        % it's a param for the estimator, hence its destined place is in params

  % TODO prob rewrite some of this stuff, at least documentation, not 
  % necessarily use currData, let them be like AFEM functions1
  
  nodeValues4e = computeNodeValues4e(s4e, u); % nodeValues4e(j) = (u(P1),u(P2),u(P3)) 
                                              %  wrt. T_j
  absNodeJumps4s = computeAbsNodeJumps4s(n4e, e4s, nodeValues4e);

  termJumps = 1/4*(...
              length4s(s4e(:, 1)).*sum(absNodeJumps4s(s4e(:, 1), :), 2)...
              + length4s(s4e(:, 2)).*sum(absNodeJumps4s(s4e(:, 2), :), 2)...
              + length4s(s4e(:, 3)).*sum(absNodeJumps4s(s4e(:, 3), :), 2)); 

  %eta4e = area4e.^(2/n).*(termF-2*alpha*(...
  %  u(s4e(:,1))'*temp1+u(s4e(:,2))'*temp2+u(s4e(:,3))'*temp3)...
  %  + alpha^2/3*area4e.*sum(u(s4e),2))... % Volumenanteil
  %  + 1/4*area4e.^(delta/n).*(...
  %   length4s(s4e(:,1))'*sum(absNodeJumps4s(s4e(:,1)),2)...
  %  +length4s(s4e(:,2))'*sum(absNodeJumps4s(s4e(:,2)),2)...
  %  +length4s(s4e(:,3))'*sum(absNodeJumps4s(s4e(:,3)),2)); %Sprunganteil


        %TODO Tien has a little error, he also has n = 2, alpha = beta = 1
        %(in documentaion of his)
  mu4e = area4e.^(2/n).*normOfDifference4e;
  % TODO sth with the volume term is pretty wrong
  xi4e = area4e.^(beta4Estimate/n).*termJumps;
  eta4e = xi4e + mu4e;
end
