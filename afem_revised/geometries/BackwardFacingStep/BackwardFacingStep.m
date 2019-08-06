function BackwardFacingStep()
% BACKWARDFACINGSTEP The backward facing step geometry with pure Dirichlet
% boundary. In the following scheme 'D' depicts Dirichlet boundary whereas
% 'N' depicts Neumann boundary. '0' denotes the point of origin.
%
%   D D D D D D D D D D D D D D D D D D D D D
%   D                                       D
%   D D D D 0                               D
%           D                               D
%           D D D D D D D D D D D D D D D D D
%
% Example:
% [c4n n4e n4sDb n4sNb] = loadGeometry('BackwardFacingStep',1);
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
