function [c4n,n4e,n4sDb,n4sNb] = computeGeometryBigSquare(red)
  c4n = [-1 -1; 1 -1; 1 1; -1 1];
  n4e = [2 4 1; 4 2 3];
  n4sDb = [4 1; 2 3; 1 2; 3 4];
  n4sNb = [];

  for k=1:red
    [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
  end
end
