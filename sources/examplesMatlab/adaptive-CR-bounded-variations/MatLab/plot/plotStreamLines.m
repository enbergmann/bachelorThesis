function plotStreamLines(c4n,n4e,u4e,startX,startY,color)
%% plotStreamLines
% plotStreamLines(c4n,n4e,u4e,startX,startY) plots streamlines
% with MATLAB function streamline and startpoints (startX,startY).
%
% (C) 2009--2013 C. Merdon, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

if nargin < 6
    color = 'k';
end
coordX = reshape(c4n(n4e,1),[],3);
coordY = reshape(c4n(n4e,2),[],3);
uX = u4e(:,1:3);
uY = u4e(:,4:6);
ax = min(c4n(:,1));
bx = max(c4n(:,1));
ay = min(c4n(:,2));
by = max(c4n(:,2));
anz = 250;
[XI,YI] = meshgrid(ax:(bx-ax)/anz:bx,ay:(by-ay)/anz:by);
UI = griddata(coordX(:),coordY(:),uX(:),XI,YI);
VI = griddata(coordX(:),coordY(:),uY(:),XI,YI);
%h = streamslice(XI,YI,UI,VI,'method','cubic');
h = streamline(XI,YI,UI,VI,startX,startY,[0.05,10000]);
set(h,'color',color);
h = streamline(XI,YI,-UI,-VI,startX,startY,[0.05,10000]);
set(h,'color',color);
end