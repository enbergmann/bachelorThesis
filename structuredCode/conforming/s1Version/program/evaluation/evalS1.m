function evalCR(path)

  if nargin < 1
    hString = 'fromServerComparison';
    path = sprintf('../../../../results/%s/conforming/zeroInitial/',...
      hString);
  end

  addpath(genpath('../'),genpath('../../../utils/'));

  dirList = strsplit(ls(path));
  dirList = dirList(not(strcmp(dirList,''))); %remove empty entries
  nDirs = length(dirList);

  nDofVec = [];
  l2ErrorVec = [];
  % l2EnergyErrorVec = [];
  
  for j=1:nDirs
    myVars = {'u','c4n','n4e','alpha','delta','red'};
    load(sprintf('%s%s/workspace.mat',path,dirList{j}),myVars{:});

  
    [~,~,n4sDb,~] = computeGeometryPolygon(red);

    nDofVec(end+1) = length(setdiff(1:size(c4n,1),unique(n4sDb)));

    uExact = @(x)gUexact(x,alpha,delta);
    l2ErrorVec(end+1) = sqrt(sum(error4eP1L2(c4n,n4e,uExact,u)));

    % gradUexact = @(x)gradG    uExact(x,delta);
    % l2EnergyErrorVec(end+1) = sqrt(sum(error4eCREnergy(c4n,n4e,gradUexact,u)));
  end

  [nDofVec,I] = sort(nDofVec);
  l2ErrorVec = l2ErrorVec(I);
  % l2EnergyErrorVec = l2EnergyErrorVec(I);
   
  name = sprintf('%s../../%sS1.txt',path,hString);
  file = fopen(name,'w');
  fprintf(file,'%s %s\n','nDof','l2Error');
  fprintf(file,'%d %.8e\n',[nDofVec;l2ErrorVec]);
  fclose(file);

  l2ErrorFig = figure('visible','off'); 
  loglog(nDofVec,l2ErrorVec);
  xlabel('nDof');
  ylabel('exact L2 error');
  name = sprintf('%s../../%sS1.png',path,hString);
  saveas(l2ErrorFig,name);

  % energyErrorFig = figure('visible','on'); 
  % loglog(nDofVec,l2EnergyErrorVec);
  % xlabel('nDof');
  % ylabel('exact energy error');
  % saveas(energyErrorFig,'energyTest1e-3.png');
 end
