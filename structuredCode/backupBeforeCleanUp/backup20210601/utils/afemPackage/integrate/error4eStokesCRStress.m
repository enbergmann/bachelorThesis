function error4e = error4eStokesCRStress(c4n,n4e,component,stressExact,uApprox,pApprox,gradUApprox)
%% error4eStokesCRStress - exact stress error of CR Stokes solution
% Compute the exact error of one row of the stress tensor for a given CR
% finite element solution of the Stokes problem.
%
% Input:    c4n         coordinates for the nodes of the mesh
%           n4e         nodes for the elements of the mesh
%           component   number of row of the stress tensor for which the
%                       error should be computed
%           stressExact exact stress function
%           uApprox     basis coefficients of the numerical solution of the
%                       velocity field w.r.t. the Crouzeix-Raviart basis
%           pApprox     basis coefficients of the numerical solution of the 
%                       pressure w.r.t. the P0 basis
%           gradUApprox corresponding row of the piecewise gradient of the
%                       discrete solution u (OPTIONAL)
%
% Output:   error4e     exact stress error on each element

    %% Initialization
    nrElems = size(n4e,1);
    s4e = computeS4e(n4e);
    
    %% Computation of discrete gradient (optional)
    if nargin < 7
        gradUApprox = zeros(nrElems,2);
        for j = 1:nrElems;
            grads = [c4n(n4e(j,:),:)'; 1 1 1] \ [-2 0; 0 -2; 0 0];
            gradUApprox(j,:) = uApprox(s4e(j,:))' * grads([3 1 2],:);
        end
    end
    
    %% Computation of discrete stress
    if component == 1
        stressUApprox = gradUApprox-[pApprox,zeros(size(gradUApprox,1),1)];
    elseif component == 2
        stressUApprox = gradUApprox-[zeros(size(gradUApprox,1),1),pApprox];
    end
    
    %% Computation of error
    error4e = sum(integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                            (stressExact(Gpts4p)-stressUApprox).^2, 6),2);
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009-2015
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
