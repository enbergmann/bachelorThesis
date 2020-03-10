function vCR = interpolationCR(currData, v)
% Computes the Crouzeix-Raviart interpolation vNC (vanishing in the midpoints
% of boundary edges) with respect to the triangulation given by [c4n, n4e, n4s]
% of the function v.
%
% interpolationCR.m
% input:  currData - struct with fields:
%                           c4n: coordinates for nodes
%                           n4s: nodes for sides
%                      length4s: lengths for sides
%         v        - 'function_handle' of the function to be interpolated
%
% output: vCR      - '(number of sides x 1)-dimensional double array' where the
%                    j-th row contains the coefficient of the CR interpolation
%                    of v wrt. the j-th side of the triangulation

  % extract necessary data
  c4n = currData.c4n;
  n4s = currData.n4s;
  length4s = currData.length4s;
  degree4Integrate = currData.degree4Integrate;

  % compute vCR
  vCR = integrate(@(n4p,Gpts4p,Gpts4ref) v(Gpts4p), ...
    c4n, n4s, degree4Integrate, length4s)./length4s;
end
