function vCR = courant2CR(n4e, v)
%% DOC 
% Computes the coefficients as a CR function of the Courant function 
% v with respect to a mesh defined by n4e
%
% courant2CR.m
% input: n4e   - nodes for elements
%        v     - '(nrNodes x 1)-dimensional double array'
%                where the j-th row contains the coefficient of the Courant
%                function w.r.t. the j-th node
%
% output: vNew - '(nrSides of the refined mesh x 1)-dimensional double array'
%                where the j-th row contains the coefficient of the CR
%                prolongation of v w.r.t. the j-th side of the new
%                triangulation

%% INIT
  s4e = computeS4e(n4e); 

%% MAIN
  % compute CR coefficients of v
  vCR = zeros(max(max(s4e)), 1);
  for elem = 1:size(n4e, 1)
    temp = v(n4e(elem, :)); % values in nodes of current element
    vCR(s4e(elem, :)) = (temp + temp([2 3 1]))/2; 
     % average the values in nodes of the current edge for the value in
     % midpoint 
     % no addition needed since CR continuous in the midpoints BUT this
     % overrides already known values for every inner edge with the same value,
     % which is fine since the Courant function v is globally continuous
  end
end
