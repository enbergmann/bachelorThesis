function plotP2(c4n, n4e, x, OPTtitle, OPTnRef)
%% Draw a P2-function.
%   plotP2(c4n, n4e, x, OPTtitle) draws the P2-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x).
%                                 The input argument OPTtitle is
%                                 optional, it sets the title of the figure.
%                                 The default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

  if nargin >= 5
    nRef = OPTnRef;
  else
    nRef = 1;
  end

  %% INITIALIZATION
  nComp = size(x, 2);
  nElem = size(n4e, 1);
  nNodes = size(c4n, 1);
  s4e = computeS4e(n4e);

  % prolongation to finer mesh
  x4e = reshape(x([n4e(:); nNodes+s4e(:)],:), nElem, 6*nComp);
  [c4n, n4e, ~, ~, x4e] =...
        P2InterpolationOnRedRefinement(c4n, n4e, zeros(0,2), zeros(0,2),...
                                       x4e, nRef);
  s4e = computeS4e(n4e);
  [~, Ind] = unique([n4e, size(c4n, 1)+s4e]);
  x = zeros(size(c4n, 1)+max(s4e(:)), nComp);
  for j = 1:nComp
    aux = x4e(:,1+6*(j-1):6*j);
    x(:,j) = aux(Ind);
  end

  %% PLOT
  if nComp == 1
    trisurf(n4e, c4n(:,1), c4n(:,2), x(1:size(c4n, 1)), 'EdgeColor', 'none');
  elseif nComp == 2
    set(gca, 'XLim', 1.1*[min(c4n(:,1)), max(c4n(:,1))],...
             'YLim', 1.1*[min(c4n(:,2)), max(c4n(:,2))]);
    n4s = computeN4s(n4e);
    mid4s = computeMid4s(c4n, n4s);
    quiver2([c4n(:,1); mid4s(:,1)], [c4n(:,2); mid4s(:,2)],...
            x(:,1), x(:,2), 'n=', 0.1, 'w=', [1 1]);
  else
    error('Invalid shape of coefficient vector');
  end

  if nargin >= 4
    title(OPTtitle);
  else
    title('');
  end

  drawnow;
end


function [c4n, n4e, n4sDb, n4sNb, x] =...
          P2InterpolationOnRedRefinement(c4n, n4e, n4sDb, n4sNb, x, iter)
% compute P2 interpolation on red-refinement of current triangulation of
% given P2 function on coarser triangulation

%% INITIALIZATION
nComp = fix(size(x, 2) / 6);

%% CONSTANTS
% TRANSFORMATION MATRICES
% child 1
T1=[1,      0,      0,     0,    0,    0;
    0,      0,      0,     0,    0,    1;
    0,      0,      0,     0,    1,    0;
    0,     -0.125, -0.125, 0.25, 0.5,  0.5;
    0.375,  0,     -0.125, 0,    0.75, 0;
    0.375, -0.125,  0,     0,    0,    0.75];
% child 2
T2=[ 0,     0,      0,     0,    0,    1;
     0,     1,      0,     0,    0,    0;
     0,     0,      0,     1,    0,    0;
     0,     0.375, -0.125, 0.75, 0,    0;
    -0.125, 0,     -0.125, 0.5,  0.25, 0.5;
    -0.125, 0.375,  0,     0,    0,    0.75];
% child 3
T3=[ 0,      0,     0,     0,    1,    0;
     0,      0,     0,     1,    0,    0;
     0,      0,     1,     0,    0,    0;
     0,     -0.125, 0.375, 0.75, 0,    0;
    -0.125,  0,     0.375, 0,    0.75, 0;
    -0.125, -0.125, 0,     0.5,  0.5,  0.25];
% child 4
T4=[ 0,      0,      0,     1,    0,    0;
     0,      0,      0,     0,    1,    0;
     0,      0,      0,     0,    0,    1;
     0,     -0.125, -0.125, 0.25,  0.5,  0.5;
    -0.125,  0,     -0.125, 0.5,  0.25, 0.5;
    -0.125, -0.125,  0,     0.5, 0.5,  0.25];


for count = 1:iter
  %% RED-REFINEMENT
  [c4n, n4e, n4sDb, n4sNb, ~, ~, ~, ~, childnr] = ...
    refineUniformRed(c4n, n4e, n4sDb, n4sNb);

  %% COMPUTATION
  x_fine = zeros(size(n4e, 1), 6*nComp);
  for j = 1:nComp
    x_fine(childnr==1,1+6*(j-1):6*j) = x(:,1+6*(j-1):6*j)*T1';
    x_fine(childnr==2,1+6*(j-1):6*j) = x(:,1+6*(j-1):6*j)*T2';
    x_fine(childnr==3,1+6*(j-1):6*j) = x(:,1+6*(j-1):6*j)*T3';
    x_fine(childnr==4,1+6*(j-1):6*j) = x(:,1+6*(j-1):6*j)*T4';
  end

  x = x_fine;
end
end
