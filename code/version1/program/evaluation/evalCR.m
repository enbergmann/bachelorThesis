function evalCR(path)

  if nargin < 1
    path = '../../../../results/forPlot_h=1e-6/nonconforming/zeroInitial/';
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
    uExact = @(x)gUexact(x,alpha,delta);
    [~,~,n4sDb,n4sNb] = computeGeometryPolygon(red);
    l2ErrorVec(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));
    nrSides = length(computeN4s(n4e));
    nDofVec(end+1) = length(computeDof(n4e,nrSides,n4sDb,n4sNb));
  end

  [nDofVec,I] = sort(nDofVec);
  l2ErrorVec = l2ErrorVec(I);


  % further plots  
  errorL2Fig = figure('visible','on'); 
  loglog(nDofVec,l2ErrorVec);
  xlabel('nDof');
  ylabel('exact L2 error');
  saveas(errorL2Fig,'test1e-6.png');
 end
