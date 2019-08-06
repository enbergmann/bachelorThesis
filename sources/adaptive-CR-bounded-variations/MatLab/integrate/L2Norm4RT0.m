function [val,valSq4e]=L2Norm4RT0(c4n,n4e,sigma4e) 
%% L2Norm4RT0
% computes the L2 norm of a RT0 function in arbitrary components
%
% part of AFEM package
% by Philipp Bringmann

  %% INITIALIZATION
  nrElem=size(n4e,1);
  nrComp=size(sigma4e,2)/3;
  area4e=computeArea4e(c4n,n4e); 

  % allocation of memory
  valSq4e=zeros(nrElem,1);

  %% COMBINATION OF COORDINATES
  % c4nT=c4n'; c4e=reshape(c4nT(:,n4e'),6,nrElem)';
  c4e=[reshape(c4n(n4e,1),nrElem,3),reshape(c4n(n4e,2),nrElem,3)];
  recombined=sum(c4e.*c4e(:,[2,3,1,5,6,4]),2);

  %% INTEGRATION
  % formula by representation with barycentric coordinates
  for k=0:nrComp-1
    valSq4e=valSq4e + (sigma4e(:,1+k*3).^2+sigma4e(:,2+k*3).^2) + 1/6*sigma4e(:,3+k*3).*...
            ( 4*sigma4e(:,1+k*3).*sum(c4e(:,1:3),2) + 4*sigma4e(:,2+k*3).*sum(c4e(:,4:6),2) +...
              sigma4e(:,3+k*3).*(sum(c4e.^2,2)+recombined) );
  end
  valSq4e=valSq4e.*area4e;
  val=sqrt(sum(valSq4e(:)));

end
