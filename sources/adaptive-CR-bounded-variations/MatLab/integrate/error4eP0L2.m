function [errorSigmaP0Sq4e, normP0funcL2Sq4e, normExfuncL2Sq4e] = error4eP0L2(c4n,n4e,P0func4e,Exfunc,intDegree)
  
  area4e = computeArea4e(c4n,n4e);
  nrElems = size(n4e,1);
  
  normP0funcL2Sq4e = area4e.*reshape(sum(sum((P0func4e.*P0func4e))),[nrElems,1,1]);
  normExfuncL2Sq4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                         permute(sum(sum(Exfunc(Gpts4p).^2)),[3 2 1]), intDegree^2,[nrElems,1,1],area4e);
  
  ExfuncTmp = @(x) permute(Exfunc(x),[3,1,2]);
  intExfunc4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                            ExfuncTmp(Gpts4p),...
                            intDegree+1,[nrElems,2,2],area4e);
  intExfunc4e = permute(intExfunc4e,[2,3,1]);
  
  mixedTerm4e = sum(sum(P0func4e.*intExfunc4e));
  mixedTerm4e = permute(mixedTerm4e,[3,2,1]);
  
  errorSigmaP0Sq4e = normP0funcL2Sq4e - 2*mixedTerm4e + normExfuncL2Sq4e;
end
