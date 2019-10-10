%TODO this folder is for functions that don't need to be documented
%TODO very old version, needs to be rewritten if needed again
% WON'T function for now
function saveScreenshot(firstScreenshot)
  dirName = sprintf(...
    '../../../../results/tvRegPrimalDualScreenshots/%s/%s',...
    firstScreenshot,datestr(now,'yy_mm_dd_HH_MM_SS'));
  mkdir(dirName);

  limitAitken = ...
    (energyVec(end)*energyVec(end-2)-energyVec(end-1)^2)/...
    (energyVec(end-2)-2*energyVec(end-1)+energyVec(end));

  name = sprintf('%s/AitkenLimitEnergyVec.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, 'Aitken Limit\n%.8g\n',limitAitken);
  fprintf(file, '%.8g\n',energyVec);
  fclose(file);

  % Use Aitken Extrapolation as limit
  energyRateAitkenFig = figure('visible','off');
  loglog(energyVec-limitAitken);
  ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
  fName = sprintf('%s/EnergyRateFigRed%dAitken.png',dirName,red);
  saveas(energyRateAitkenFig,fName);

  energyRateAitkenFigSemilog = figure('visible','off');
  semilogy(energyVec-limitAitken);
  ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
  fName = sprintf('%s/EnergyRateFigSemilogRed%dAitken.png',dirName,red);
  saveas(energyRateAitkenFigSemilog,fName);

  % Use last energyVec value as limit
  energyRateFig = figure('visible','off');
  loglog(energyVec-energyVec(end));
  ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
  fName = sprintf('%s/EnergyRateFigRed%d.png',dirName,red);
  saveas(energyRateFig,fName);

  energyRateFigSemilog = figure('visible','off');
  semilogy(energyVec-energyVec(end));
  ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
  fName = sprintf('%s/EnergyRateFigSemilogRed%d.png',dirName,red);
  saveas(energyRateFigSemilog,fName);
end
