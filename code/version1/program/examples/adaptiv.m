function adaptiv(red,terminate)

  addpath(genpath('../'),genpath('../../../utils/'));

  figVisible = 'on';
  % set(0,'DefaultFigureVisible','off');
  
  initalU = 'zero'; %f, zero

  miscMsg = 'adaptiv algorithmn';
  expName = 'adaptivTest';
  
  tau = 1/2;
  % tau = .1;
  h = 2^(-red);

  alpha = 1; 
  delta = 1;     
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red);
  
  n4s = computeN4s(n4e);
  
  %% given analytic example

  f=@(x)g(x,alpha,delta);  
  % uExact=@(x)gUexact(x,alpha,delta);  
  
  %% f = 0
  %f=@(x)0;  
  %uExact=@(x)0;  
  
  %mid4s = computeMid4s(c4n,n4s);
  %

  if strcmp(initalU,'zero')
    u = zeros(size(n4s,1),1);
    message = sprintf('on unit circle, inital u = 0');
    dirInfoName = sprintf('zeroInitial');
  elseif strcmp(initalU,'f')
    u = interpolationNC(f,c4n,n4e,n4s);
    message = sprintf('on unit circle, inital u = I_NC(f)');
    dirInfoName = sprintf('fInitial');
  end
  
  du = computeGradientNC(c4n,n4e,u);
  Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
  Lambda(isinf(Lambda)) = 0;
  Lambda(isnan(Lambda)) = 0;
  
  %  Lambda = zeros(size(n4e,1),2);
  %  u = zeros(size(n4s,1),1);
  
  %% MAIN AFEM LOOP

  minNrDoF = 1000;
  eta4nrDoF = sparse(1,1);

  while( true )
    % SOLVE
    tic;
    [u,corrVec,energyVec,nrDof] = ...
      tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u,Lambda);
    time = toc; 
      % [x,nrDoF] = solveCRPoisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
   
    saveResults('CR',expName,dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
      terminate,time,corrVec,energyVec,tau,miscMsg);

    %TODO
    % ESTIMATE
    [eta4s,n4s] = estimateCREtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
    eta4nrDoF(nrDoF) = sqrt(sum(eta4s));
    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
        '; estimator = ',num2str(eta4nrDoF(nrDoF))]);

    % Check Termination
    if nrDoF >= minNrDoF, break, end;

    % MARK
    n4sMarked = markBulk(n4s,eta4s);

    % REFINE
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);

    %% Plot mesh, solution and convergence graph.
    % TODO add to save results things like plotTriangulation and so on
    % and stop saving workspace elemnts that are too big (c4n), update
    % nrDof to be sth returned so that it doesn't have to be computed
    % in evaluate
    figure;
    plotTriangulation(c4n,n4e);
    figure;
    plotCR(c4n,n4e,x,{'CR Solution'; [num2str(nrDoF) ' degrees of freedom']});
    nrDoF4lvl = find(eta4nrDoF);
    eta4lvl = eta4nrDoF(nrDoF4lvl);
    figure;
    plotConvergence(nrDoF4lvl,eta4lvl,'\eta_l');
  end
end
