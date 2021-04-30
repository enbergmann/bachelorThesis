function l1NormOfJump4s = computeL1NormOfJump4s(currData, output)
%% DOC
% Computes for a CR function u the L^1 norm of the jump ||[u]_F||_{L^1(F)} for
% all F\in\Fcal w.r.t. a triangulation given by [c4n, n4e].
%
% computeL1NormOfJump4s.m
% input: currData - 'struct' with fields:
%                          n4e: nodes for elements 
%                          n4s: nodes for sides
%                      nrSides: number of sides
%                          s4e: sides for elements
%                     length4s: lengths for sides
%                          e4s: elements for sides
%
%        output   - 'struct' with fields: 
%                     u: '(nrSides x 1)-dimensional double array' where the
%                        j-th row contains the coefficient of a CR function u
%                        w.r.t. the j-th edge of the triangulation 
%
% output: l1NormOfJump4s - '(nrSides x 1)-dimensional double array' where the
%                          j-th row contains the L1 norm of the jump of u across
%                          the j-th side of the triangulation (if a side is an
%                          outer edge the jump is defined as the restriction
%                          of u to this side)

%% INIT
  % extract necessary information from currData
  n4e = currData.n4e;
  n4s = currData.n4s;
  nrSides = currData.nrSides;
  s4e = currData.s4e;
  length4s = currData.length4s;
  e4s = currData.e4s;

  % extract necessary information from output
  u = output.u;


  % get the absolute jumps in the two nodes for of each side in the
  % triangulation
  nodeValues4e = computeNodeValuesCR4e(s4e, u); 
  absNodeJumps4s = computeAbsNodeJumps4s(n4e, n4s, e4s, nrSides, nodeValues4e);

%% MAIN
  l1NormOfJump4s = 1/4*length4s.*sum(absNodeJumps4s, 2);
    % cf. thesis
end
