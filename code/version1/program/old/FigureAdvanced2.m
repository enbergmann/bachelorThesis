function Fig = FigureAdvanced2(OldFig)
if nargin == 1
%% Replace existing Figure
% clear contents
clf(OldFig);
% make current figure
figure(OldFig);
% rename for output
Fig = OldFig;
else
%% Create new figure
FigCount = length(findobj('type','figure'));
Fig = figure(FigCount + 1);
%replace new figure
ScreenInfo = get(0,'ScreenSize');
if ~isequal(ScreenInfo,[1 1 1 1])
ScreenWidth = ScreenInfo(3);
ScreenHeight = ScreenInfo(4);

FigWidth = 560;
FigHeiht = 420;
Border = 0;

FiguresPerRow = floor((ScreenWidth - Border)/(FigWidth+Border));
FiguresPerColumn = floor((ScreenHeight - Border)/(FigHeiht+Border));

FigRow = ceil(Fig/FiguresPerRow);
FigColumn = Fig - (FigRow-1)*FiguresPerRow;

set(Fig,'Position',[FigColumn*Border + (FigColumn - 1)*FigWidth,...
ScreenHeight - FigRow*Border - FigRow*(FigHeiht),...
FigWidth, FigHeiht]);
end
end
end