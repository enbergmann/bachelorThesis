function ENew = computeEnergy(area4e,u,du,alpha,temp,MAMANC)
  ENew =  area4e'*sqrt(sum(du.^2,2)) ... % L^1 norm uf gradNC u
               - u'*temp ... % integral f u
               + alpha/2* u'*MAMANC*u;    % L^2 norm uf u
end
