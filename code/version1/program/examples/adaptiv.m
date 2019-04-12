function adaptiv(red,terminate)

  addpath(genpath('../'),genpath('../../../utils/'));

  figVisible = 'off';
  % set(0,'DefaultFigureVisible','off');
  
  initalU = 'zero'; %f, zero

  if strcmp(initalU,'zero')
    message = sprintf('on unit circle, inital u = 0');
    dirInfoName = sprintf('zeroInitial/%s',datestr(now,'yy_mm_dd_HH_MM_SS'));
  elseif strcmp(initalU,'f')
    message = sprintf('on unit circle, inital u = I_NC(f)');
    dirInfoName = sprintf('fInitial',datestr(now,'yy_mm_dd_HH_MM_SS'));
  else
    message = sprintf('on unit circle, unknownInitialF');
    dirInfoName = sprintf('unknownInitialF',datestr(now,'yy_mm_dd_HH_MM_SS'));
  end

  miscMsg = 'adaptiv algorithmn';
  expName = 'adaptivTest';
  
  tau = 1/2;
  % tau = .1;
  h = 2^(-red);

  alpha = 1; 
  delta = 1;     
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red);
  polygonGeometry = true;
  
  
  %% given analytic example

  % f=@(x)0;  
  f=@(x)g(x,alpha,delta);  
  uExact=@(x)gUexact(x,alpha,delta);  
  
  %% f = 0
  %f=@(x)0;  
  %uExact=@(x)0;  
  
  %mid4s = computeMid4s(c4n,n4s);
  %
  
  %% MAIN AFEM LOOP

  minNrDoF = 10000000;
  eta4lvl = [];
  nrDoF4lvl = [];
  l2Error4lvl = [];

  while( true )
    % SOLVE
    
    n4s = computeN4s(n4e);

    if strcmp(initalU,'zero')
      u = zeros(size(n4s,1),1);
    elseif strcmp(initalU,'f')
      u = interpolationNC(f,c4n,n4e,n4s);
    end
    
    du = computeGradientNC(c4n,n4e,u);
    Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
    Lambda(isinf(Lambda)) = 0;
    Lambda(isnan(Lambda)) = 0;
    
    %  Lambda = zeros(size(n4e,1),2);
    %  u = zeros(size(n4s,1),1); 
    tic;
    [u,corrVec,energyVec,nrDof] = ...
      tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u,Lambda);
    time = toc; 
    nrDoF4lvl(end+1) = nrDof;
   
    % ESTIMATE
    eta4e = estimateError4e(u,f,c4n,n4e,n4sDb,n4sNb,alpha,delta)
    eta4lvl(end+1) = sum(eta4e);
    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDof),...
      '; estimator = ',num2str(eta4lvl(end))]); 

    l2Error4lvl(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));

    saveResults('CR',expName,dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
      terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
      true,nrDoF4lvl,eta4lvl,l2Error4lvl);

    % Check Termination
    if nrDof >= minNrDoF, break, end;

    % MARK
    n4sMarked = markBulk(n4e,eta4e,0.5);

    % REFINE
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);

    if polygonGeometry
      temp=unique(n4sDb);
      c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
    end
  end
end
