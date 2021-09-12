function vGradCR = gradientCR(currData, v)
%% DOC
% Computes the piecewise gradient of the Crouzeix-Raviart function v (vanishing
% in the midpoints of boundary edges) with respect to the triangulation given
% by [c4n, n4e].
%
% gradientCR.m
% input: currData - 'struct' with fields:
%                       nrElems: number of elements
%                           s4e: sides for elements
%                     gradsCR4e: '(3 x 2 x nrElems)-dimensional double array'
%                                where the j-th row of the k-th matrix contains
%                                the gradient of the local CR-basis function
%                                w.r.t. the j-th edge of the k-th element
%        v        - '(nrSides x 1)-dimensional double array' where the j-th
%                   row contains the coefficient of v w.r.t. the j-th side of
%                   the triangulation 
%
% output: vGradCR - '(nrElems x 2)-dimensional double array' where the j-th
%                   row contains the gradient of v on the j-th triangle 

%% INIT
  % extract necessary information from currData
  nrElems = currData.nrElems;
  s4e = currData.s4e;
  gradsCR4e = currData.gradsCR4e;

%% MAIN
  vReshaped = reshape(v(s4e)', 3*nrElems, 1);
  gradsReshaped = reshape(permute(gradsCR4e, [1 3 2]), nrElems*3, 2);
    % for the CR coeffiecient v(k,j) and the gradient (as a row) g(k,j) of the
    % local CR basis function w.r.t. the j-th edge of the k-th triangle of the
    % triangulation it holds
    %   vReshaped = [v(1,1); v(1,2); v(1,3); ...
    %                v(2,1); v(2,2); v(2,3); ...
    %                ...
    %                v(nrElems,1); v(nrElems,2); v(nrElems,3)] and
    %   gradsReshaped = [g(1,1); g(1,2); g(1,3); ...
    %                    g(2,1); g(2,2); g(2,3); ...
    %                    ...
    %                    g(nrElems,1); g(nrElems,2); g(nrElems,3)]
  gradsWeighted = vReshaped.*gradsReshaped;
    %   gradsWeighted = ...
    %     [v(1,1)*g(1,1); v(1,2)*g(1,2); v(1,3)*g(1,3); ...
    %      v(2,1)*g(2,1); v(2,2)*g(2,2); v(2,3)*g(2,3); ...
    %      ...
    %      v(nrElems,1)*g(nrElems,1); ...; v(nrElems,3)*g(nrElems,3)]
  vGradCR = reshape(sum(reshape(gradsWeighted, 3, nrElems*2)), nrElems, 2);
    % add every three rows of gradsWeighted to obtain
    %   vGradCR(k) = [v(k,1) v(k,2) v(k,3)]*[g(k, 1); g(k, 2); g(k, 3)]
end
