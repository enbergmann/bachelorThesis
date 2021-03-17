function [errorDivU4e,normDivUSq4e,normDivUexactSq4e] = error4eDivP1(geometryData,localData,solutionData,exSolutionData,intDegree)
   area4e    = localData.area4e; %localData{15};
   divU4e    = solutionData.divU4e; %solutionData{8};
   divUexact = exSolutionData.divUexact; %exSolutionData{10};
   c4n       = geometryData.c4n; %geometryData{1};
   n4e       = geometryData.n4e; %geometryData{2};
   nrElems   = size(n4e,1);
   
   normDivUSq4e = area4e.*divU4e.^2;
   
   mixedTerm4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref) divUexact(Gpts4p),...
                                     intDegree,[nrElems,1]);
   
   mixedTerm4e = divU4e.*mixedTerm4e;
   
   normDivUexactSq4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref) divUexact(Gpts4p).^2,...
                                     2*intDegree,[nrElems,1]);
   
   errorDivU4e = normDivUSq4e - 2*mixedTerm4e + normDivUexactSq4e;
end
