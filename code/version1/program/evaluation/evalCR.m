function evalCR(path)

  if nargin < 1
    path = '../../../../results/forPlot_h=1e-3/nonconforming/zeroInitial/';
  end

  addpath(genpath('../'),genpath('../../../utils/'));

  dirList = strsplit(ls(path));
  dirList = dirList(not(strcmp(dirList,''))); %remove empty entries
  nDirs = length(dirList);

  nDofVec = [];
  l2ErrorVec = [];
  l2EnergyErrorVec = [];
  
  for j=1:nDirs
    myVars = {'u','c4n','n4e','alpha','delta','red'};
    load(sprintf('%s%s/workspace.mat',path,dirList{j}),myVars{:});

    [~,~,n4sDb,n4sNb] = computeGeometryPolygon(red);

    nrSides = length(computeN4s(n4e));
    nDofVec(end+1) = length(computeDof(n4e,nrSides,n4sDb,n4sNb));

    uExact = @(x)gUexact(x,alpha,delta);
    l2ErrorVec(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));

    gradUexact = @(x)gradGuExact(x,delta);
    l2EnergyErrorVec(end+1) = sqrt(sum(error4eCREnergy(c4n,n4e,gradUexact,u)));
  end

  [nDofVec,I] = sort(nDofVec);
  l2ErrorVec = l2ErrorVec(I);
  l2EnergyErrorVec = l2EnergyErrorVec(I);


  % further plots  
  errorL2Fig = figure('visible','on'); 
  loglog(nDofVec,l2ErrorVec);
  xlabel('nDof');
  ylabel('exact L2 error');
  saveas(errorL2Fig,'l2Test1e-3.png');

  energyErrorFig = figure('visible','on'); 
  loglog(nDofVec,l2EnergyErrorVec);
  xlabel('nDof');
  ylabel('exact energy error');
  saveas(energyErrorFig,'energyTest1e-3.png');
 end
