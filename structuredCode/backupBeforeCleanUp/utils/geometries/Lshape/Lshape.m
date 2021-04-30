function Lshape()
% LSHAPE The L-Shape geometry with pure Dirichlet boundary. In the 
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
