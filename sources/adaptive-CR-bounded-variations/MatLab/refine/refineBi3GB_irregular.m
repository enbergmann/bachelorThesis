function [c4nNew, n4eNew, n4sDbNew, n4sNbNew, err4eNew] = ...
                refineBi3GB_irregular(c4n, n4e, n4sDb, ...
                n4sNb, n4sMarked, err4e, rhsf, regular)
		
 %% refineBi3GB_irregular
 %		- Refine using the Bisec3-Green-Blue-strategy
 %                furthermore calculate err4eNew based on rhsf 
 %                used in approx based on err4e.
 %                	For details on data structures 
 %		 	(cf. refineBi3GB.m).
 %  However, do not use closure or completion if not regular. Thus, the 
 %  output triangulation will probably has hanging nodes.
 %
 % Input: c4n, n4e, n4sDb, n4sNb 
 %                       - triangulation 
 %        n4sMarked      - set of marked sides to be refined 
 %        err4e = [rhsf4e, etilde4e] 
 %                       - error functional and error functional 
 %                         etilde from Tresholding second algorithm
 %                         of the input triangulation
 %        rhsf           - function handle of the error 
 %                         functional for TSA 
 %       regular         - true: run closure after each marking step
 %                       - false : do not run closure simultaneaously
 %                              i.e., after each marking loop
 
 % Output: c4nNew, n4eNew, n4sDbNew, n4sNbNew  
 %                       - new triangulation
 %         err4eNew = [rhsf4eNew, etilde4eNew]
 %                       - error functional and error functional 
 %                         etilde from Tresholding second algorithm
 %                         of the output triangulation

%    Copyright (C) 2009  Carsten Carstensen; Numerical Analysis Group
%    Copyright (C) 2013  Hella Rabus
%                   
%    You should have received LICENSE.txt along with this file
%    that gives further information an the license.
%    See the GNU General Public License for more details.

 
 %% Initialisation
 nrNodes = size(c4n,1);
 nrElems = size(n4e,1);
 %%  Avoid Closure Algorithm ?
 if (regular)
    n4sRefine = closure(n4e,n4sMarked);
 else
    n4sRefine = n4sMarked;
 end
 
 %% Compute newNodes4n to find new nodes faster.
 newNodes4s = sparse(n4sRefine(:,1),n4sRefine(:,2),...
                    (1:size(n4sRefine,1))'+ nrNodes, ...
                        nrNodes, nrNodes);
 if (regular) 
     newNodes4s = newNodes4s + newNodes4s';
 end

 %% Compute coordinates of new nodes.
 mid4sRefine = (c4n(n4sRefine(:,1),:)+c4n(n4sRefine(:,2),:))/2;
 c4nNew = [c4n;mid4sRefine];

 %% bisec3 refinement
 % Count elements in new triangulation. For each newNode, one
 % new elements will be created.
 if (regular)
    nrNewElems = nrElems+2*size(n4sRefine,1)-size(n4sDb,1)...
 			-size(n4sNb,1);
 else
    nrNewElems = nrElems+nnz(newNodes4s);
 end
 n4eNew = zeros(nrNewElems,3);
 err4eNew = zeros(nrNewElems,2);
 % index to keep track of the current element number in n4eNew
 ind = 0; 
 for curElem = 1 : nrElems
     curNodes = n4e(curElem,:);
     curErr = err4e(curElem,:);
     curNewNodes = [newNodes4s(curNodes(1),curNodes(2));
                    newNodes4s(curNodes(2),curNodes(3));
                    newNodes4s(curNodes(3),curNodes(1));
                   ];
     nrNewNodes4curElem = nnz(curNewNodes);
     if nrNewNodes4curElem == 0 % no refinement
         n4eNew(ind+1,:) = curNodes;
         err4eNew(ind+1,:) = curErr;
         ind = ind+1;
     elseif nrNewNodes4curElem == 1 % green refinement
         n4eNew(ind+1:ind+2,:) = ...
             [ curNodes(3)    curNodes(1) curNewNodes(1);
               curNodes(2)    curNodes(3) curNewNodes(1);
             ];
         % Compute rhsf_term
         rhs4e = rhsf(c4nNew, n4eNew(ind+1:ind+2,:));
         % Compute etilde
         if sum(curErr)<=0
     	    err4eNew(ind+1:ind+2,:) = [rhs4e zeros(2,1)];
         else
     	    err4eNew(ind+1:ind+2,:) = [rhs4e ...
            sum(rhs4e)/sum(curErr)*curErr(2)*ones(2,1)];
         end
         ind = ind+2;
     elseif nrNewNodes4curElem == 2
         if curNewNodes(2) > 0 % blue right
             curn4eNew = ...
                 [ curNodes(3)    curNodes(1)    curNewNodes(1);
                   curNewNodes(1) curNodes(2)    curNewNodes(2);
                   curNodes(3)    curNewNodes(1) curNewNodes(2);
                   curNodes(2)    curNodes(3)    curNewNodes(1)
                 ];
             n4eNew(ind+1:ind+3,:) = curn4eNew(1:3,:);
             % Compute rhsf_term
             rhs4e = rhsf(c4nNew, curn4eNew);
             % Compute etilde
             if sum(curErr)<=0
                err4T1=0;
             else
                err4T1 = (rhs4e(1)+rhs4e(4))/...
                    sum(curErr)*curErr(2);
             end
             err4eNew(ind+1:ind+3,:) = [rhs4e(1:3) [err4T1 ;...
                 sum(rhs4e(2:3),1)/(rhs4e(4)+err4T1)*...
                 err4T1*ones(2,1)]];

         else % blue left
             curn4eNew = ...
                 [ curNodes(1)    curNewNodes(1) curNewNodes(3);
                   curNewNodes(1) curNodes(3)    curNewNodes(3);
                   curNodes(2)    curNodes(3)    curNewNodes(1);
                   curNodes(3)    curNodes(1)    curNewNodes(1)
                 ];
             n4eNew(ind+1:ind+3,:) = curn4eNew(1:3,:);
             % Compute rhsf_term
             rhs4e = rhsf(c4nNew, curn4eNew);
             % Compute etilde
             if sum(curErr)<=0
                 err4T3=0;
             else 
                 err4T3 = (rhs4e(3)+rhs4e(4))/...
		 		sum(curErr)*curErr(2);
             end
             err4eNew(ind+1:ind+3,:) = [rhs4e(1:3) ...
                 [sum(rhs4e(1:2),1)/(rhs4e(4)+err4T3)*...
                 err4T3*ones(2,1); err4T3]];
         end
         ind = ind+3;
     elseif nrNewNodes4curElem == 3 % bisec3 refinement
         curn4eNew = ...
             [ curNewNodes(1) curNodes(2)    curNewNodes(2);
               curNodes(3)    curNewNodes(1) curNewNodes(2);
               curNodes(1)    curNewNodes(1) curNewNodes(3);
               curNewNodes(1) curNodes(3)    curNewNodes(3);
               curNodes(3)    curNodes(1)    curNewNodes(1);
               curNodes(2)    curNodes(3)    curNewNodes(1)
             ];
         n4eNew(ind+1:ind+4,:) = curn4eNew(1:4,:);
         % Compute rhsf_term
         rhs4e = rhsf(c4nNew, curn4eNew);
         % Compute etilde
         if sum(curErr)<=0
             err4T56 = 0;
         else
             err4T56 = sum(rhs4e(5:6),1)/sum(curErr)*curErr(2);
         end
         err4eNew(ind+1:ind+4,:) = [rhs4e(1:4) ...
             [ sum(rhs4e(1:2),1)/(rhs4e(6)+err4T56)*...
	         err4T56*ones(2,1); ...
	         sum(rhs4e(3:4),1)/(rhs4e(5)+err4T56)*...
	         err4T56*ones(2,1)...
	     ]];
         ind = ind+4;
     end
 end

 %% refinement of  Dirichlet boundary
 nrNewDbSides= size(intersect(sort(n4sDb,2),sort(n4sRefine,2),...
 			'rows'),1);
 n4sDbNew = zeros(nrNewDbSides,2);
 ind = 0;
 for curSide = 1 : size(n4sDb,1)
     curNodes = n4sDb(curSide,:);
     curNewNodes = newNodes4s(curNodes(1),curNodes(2));
     if curNewNodes == 0
         n4sDbNew(ind + 1,:) = curNodes;
         ind = ind + 1;
     else
         n4sDbNew(ind + 1:ind + 2,:) = ...
             [ curNodes(1)    curNewNodes;
               curNewNodes    curNodes(2);
             ];
         ind = ind + 2;
     end
 end

 %% refinement of  Neumann boundary
 nrNewNbSides= size(intersect(sort(n4sNb,2),sort(n4sRefine,2),...
 			'rows'),1);
 n4sNbNew = zeros(nrNewNbSides,2);
 ind = 0;
 for curSide = 1 : size(n4sNb,1)
     curNodes = n4sNb(curSide,:);
     curNewNodes = newNodes4s(curNodes(1),curNodes(2));
     if curNewNodes == 0
         n4sNbNew(ind + 1,:) = curNodes;
         ind = ind + 1;
     else
         n4sNbNew(ind + 1:ind + 2,:) = ...
             [ curNodes(1)    curNewNodes;
               curNewNodes    curNodes(2);
             ];
         ind = ind + 2;
     end
 end
 if sum(isnan(err4eNew))>0
     error('NaN tritt in err4eNew auf');
 end
end
