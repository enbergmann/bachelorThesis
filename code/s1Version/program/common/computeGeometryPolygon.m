function [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(j)

c4n = [0 0; 0 1; -1 0; 0 -1; 1 0];
n4e=[1 2 3; 1 3 4; 1 4 5; 1 5 2];
n4sDb=n4e(:,[2,3]);
n4sNb=[];

for k=1:j
    [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
    temp=unique(n4sDb);
    c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
end
%figure; plotTriangulation ( c4n, n4e );
n4sNb=[];
end
