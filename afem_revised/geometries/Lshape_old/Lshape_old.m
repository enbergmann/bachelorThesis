function Lshape_old()
% LSHAPE_OLD The L-Shape geometry with pure Dirichlet boundary. In the
% following scheme 'D' depicts Dirichlet boundary whereas 'N' depicts
% Neumann boundary.
%
%     DDDDDDDDDDDDD
%     D           D
%     D           D
%     D     DDDDDDD
%     D     D
%     D     D
%     DDDDDDD
%
% Example:
% [c4n n4e n4sDb n4sNb] = loadGeometry('Lshape',1);
%
% See also LSHAPENB, LSHAPEROT
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
