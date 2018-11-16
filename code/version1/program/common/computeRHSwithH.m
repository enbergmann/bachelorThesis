function [b,temp] = computeRHS(c4n,n4e,s4e,nrSides,area4e,du,tau,Lambda,...
    nrElems,temp1,temp2,temp3,h)
  %% Create the right-hand side b  (midpoint rule might be seperated from integrate, integrate wo loop)
  b = zeros(nrSides,1);
  temp = zeros(nrSides,1);
  for elem = 1 : nrElems
      gradsNC = [ones(1,3); c4n(n4e(elem,[3 1 2]),:)']\[zeros(1,2); -2*eye(2)]; % gradients for CR basis
       bLocal = (h*du(elem,:)/tau - Lambda(elem,:) ) * gradsNC';  
       %TODO here is h
      temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
      b(s4e(elem,:)) = b(s4e(elem,:)) + area4e(elem)*bLocal'; % right-hand side
      % midpoint rule
%       mid = mid4e(elem,:);     % midpoint of this element
%       b(sides) = b(sides) + area*f(mid)*ones(3,1)/3 +area*bLocal'; % right-hand side
  end
  b = b + temp;
end
