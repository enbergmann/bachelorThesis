function vCR = interpolationCR(params, currData, v)
%% DOC
% Computes the Crouzeix-Raviart interpolation vCR (vanishing in the midpoints
% of boundary edges) with respect to the triangulation given by [c4n, n4e, n4s]
% of the function v.
%
% interpolationCR.m
% input: params   - 'struct' with fields:
%                     degree4Integrate: 'uint64' up to which the integration in
%                                       integrate must be exact
%        currData - 'struct' with fields:
%                          c4n: coordinates for nodes
%                          n4s: nodes for sides
%                     length4s: lengths for sides
%        v        - 'function_handle' of the function to be interpolated
%
% output: vCR - '(nrSides x 1)-dimensional double array' where the j-th row
%               contains the coefficient of the CR interpolation of v w.r.t.
%               the j-th side of the triangulation

%% INIT
  % extract necessary parameters from params
  degree4Integrate = params.degree4Integrate;

  % extract necessary information from currData
  c4n = currData.c4n;
  n4s = currData.n4s;
  length4s = currData.length4s;

%% MAIN
  vCR = integrate(c4n, n4s, @(n4p, Gpts4p, Gpts4ref) v(Gpts4p), ...
    degree4Integrate, length4s)./length4s;
end
