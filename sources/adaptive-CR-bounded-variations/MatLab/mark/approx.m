function [c4n, n4e, n4sDb, n4sNb, err4e,c4n_nc, n4e_nc, ...
		n4sDb_nc, n4sNb_nc, err4e_nc] = ...
		approx(rhsf, c4n_nc, ...
		n4e_nc, n4sDb_nc, n4sNb_nc, err4e_nc, ...
		c4n_initial, n4e_initial, n4sDb_initial,...
		n4sNb_initial, err4e_initial, ...
		threshold)

 %% approx   - generates an optimal triangulation wrt.
 %		the reduction of rhsf, s.t. 
 %		sum(err4e(:,1)) <= threshold.
 %	     - this algorithm marks elements as in approx
 %	      (cf. BBdV); completion is applied after tsa has 
 %	      finished.
 % INPUT                 
 % 	rhsf - error function that has to be reduced 
 %              (e.g. osc(f,T) or vol(f,T))
 %  c4n_nc, n4e_nc, n4sDb_nc, n4sNb_nc 
 %           - last output of approx but without running completion
 %	       (if Case B has been applied at least once before) 
 %             otherwise initial triangulation
 %  err4e_nc - matrix of error functional for the _nc 
 %             triangulation
 %             size(err4e_nc)=[size(n4e_nc,1),2]
 %             err4e_nc(:,1)= error functional for elements
 %	        err4e_nc(:,2)= tilde e, weighted error functional as  
 %              defined by Binev, Dahmen and deVore
 %  c4n_initial, n4e_initial, n4sDb_initial, n4sNb_initial
 %            - initial triangulation
 %  err4e_initial
 %	      - matrix of error functional of the _initial 
 %		  triangulation
 %  intDegree 
 %	      - degree for integrating when computing the 
 %		    error functional
 %  threshold - tolerated threshold for the error functional
 %
 %   OUTPUT
 %     c4n_nc, n4e_nc, n4sDb_nc, n4sNb_nc 
 %	      - refined triangulation (no completion,  
 %	        not regular) with sum(err4e(:,1)) <= threshold
 %	      - for later use, if TSA needs to be started again
 % 	c4n, n4e, n4sDb, n4sNb
 %	      - refined triangulation (after completion, 
 %	        regular) with sum(err4e(:,1)) <= threshold
 %	err4e - matrix of the error functional


%    Copyright (C) 2013  Hella Rabus
%                   
%    You should have received LICENSE.txt along with this file
%    that gives further information an the license.
%    See the GNU General Public License for more details.

 
 %% Initialisation, tolerance for marking such that all elements with
 % rhs(elem) >= max(rhs4elem)-max(rhs4elem)*rhsf_factor
 rhsf_factor = 1e-6;


 %% Computing error functional on last result of TSA 
 % (possibly not regular triangulation) 
 errT = sum(err4e_nc(:,1));
 errT_old = errT;
 
 %% Break if threshold is satisfied
 if (errT<=threshold)
     loop=false;
 else
     loop=true;
 end

 %% Reduce error functional by refining the triangulation
 %  Output + waitbar
 fprintf(1,' APPROX with threshold = %8.2e:', threshold);
 fprintf(1,', ||rhs f||^2= %8.2e', errT);
 fprintf(1,'\n   Compute opt. Triangulation: %5.2f %%, ||rhs f||^2: %8.2e', ...
 	0,errT);
while (loop)
     %% MARK elements with largest e tilde
         maxerr4e=max(err4e_nc(:,2));
         I = err4e_nc(:,2) >= (1-rhsf_factor)*maxerr4e;
     % choose the reference edge
     n4sMarked = n4e_nc(I,1:2);
     
     % markMaximum marks all edges of the triangles with maxError
     %n4sMarkedb = markMaximum(n4e_nc, err4e_nc(:,2), 1.0);
    
     %% REFINE those elements and compute [rhsf, etilde]
     [c4n_nc,n4e_nc,n4sDb_nc,n4sNb_nc, err4e_nc] = ...
         refineBi3GB_irregular(c4n_nc, n4e_nc, n4sDb_nc, ...
         n4sNb_nc, n4sMarked, err4e_nc, rhsf, false);
     errT = sum(err4e_nc(:,1));
     
     status = min((errT_old-errT)/(errT_old - threshold)*100,100);
     fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
     fprintf(1,'\b\b\b\b\b\b\b\b %6.2f%%, ||rhs f||^2: %9.3e', ...
         status, errT);

     %% BREAK if rhsf(f,T) satisfies threshold
     if errT<=threshold, loop=false; end;
 end
 fprintf(1,'\n');
 %% Completion of the triangulation to generate a regular one
 [c4n, n4e, n4sDb, n4sNb, err4e]=completion(c4n_nc,...
     c4n_initial, n4e_initial, n4sDb_initial, ...
     n4sNb_initial, err4e_initial,rhsf);
end
