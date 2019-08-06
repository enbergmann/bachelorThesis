function [c4n, n4e, n4sDb, n4sNb,n4sDb2,n4sNb2] = loadGeometry_elasticity(name, componentBoundary, OPTRefinementLevel)
%% loadGeometry - load data for a mesh.
% [c4n n4e n4sDb n4sNb] = loadGeometry('name') loads the data structures for
%   the mesh named 'name'. Optionally, the second parameter will
%   cause the mesh to be refined a given number of times using the uniform
%   red strategy.
% Example:
% [c4n n4e n4sDb n4sNb] = loadGeometry('LShape',2) loads the
%   mesh called 'LShape' and refines it two times.

    %% Load the geometry data.
    c4n = load([name,'_c4n.dat']);
    n4e = load([name,'_n4e.dat']);
    n4sDb = load([name,'_n4sDb.dat']);
    n4sNb = load([name,'_n4sNb.dat']);
    
    if componentBoundary == 1
       n4sDb2 = load([name,'_n4sDbSecond.dat']);
       n4sNb2 = load([name,'_n4sNbSecond.dat']);
    else
       n4sDb2 = load([name,'_n4sDb.dat']);
       n4sNb2 = load([name,'_n4sNb.dat']); 
    end
    
    %% Initial refinement.
    if nargin < 2
        OPTRefinementLevel = 0;
    end
    for i=1:OPTRefinementLevel
       if componentBoundary == 1
          [c4n, n4e, n4sDb, n4sNb,n4sDb2,n4sNb2] = refuneUniformRed_elasticity(c4n,n4e,n4sDb,n4sNb, n4sDb2, n4sNb2);
       else
        [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
        n4sDb2 = n4sDb;
        n4sNb2 = n4sNb;
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009 
% Numerical Analysis Group
% Prof. Dr. Carsten Carstensen
% Humboldt-University  
% Departement of Mathematics
% 10099 Berlin
% Germany
%
% This file is part of AFEM.
%
% AFEM is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% AFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
