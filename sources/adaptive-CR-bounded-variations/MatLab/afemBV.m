function afemBV(f, c4n, n4e, n4sDb, theta, minDof,...
                alpha, gamma, epsilon, tau, degreeF, identifier, uExact, DuExact)
    
    level = 1;
    if nargin > 12
        exactS = true;
    else
        exactS = false;
    end
    
    % compute exact energy
    if exactS
        c4n2 = c4n;
        n4e2 = n4e;
        n4sDb2 = n4sDb;
        for j = 1:8
            [c4n2, n4e2, n4sDb2, ~] = refineUniformRed(c4n2,n4e2,n4sDb2,zeros(0,2));
        end
        area4e2 = computeArea4e(c4n2, n4e2);
        
        % u exact is in H1
        % E(u) = alpha/2*||u||^2 + int_Omega 
        EnergyExact = ...
            alpha/2*integrate(c4n2, n4e2, @(n4p,x,x4Ref)sum(uExact(x).^2,2), 8, area4e2)...
            + integrate(c4n2, n4e2, @(n4p,x,x4Ref)sqrt(sum(DuExact(x).^2,2)), 8, area4e2)...
            - integrate(c4n2, n4e2, @(n4p,x,x4Ref)(f(x).*uExact(x)), 8, area4e2);
        EnergyExact = sum(EnergyExact);
    end
    
    folderName = ['results/', identifier, '/theta',ceil(theta*100),'Dof',num2str(minDof)];
    
    while true
        %% SOLVE ------------------------------------------------------------------------------
        % initiate primal dual algorithm
        [STIMACR, MASSCR, intFSq4e, intFCR14e, intFCR24e, intFCR34e,...
            gradsCR4e, s4e, dof, nrSides, area4e] = initPrimalDual(c4n, n4e, n4sDb, f, degreeF);
        
        nrDof = length(dof);
        e4s = computeE4s(n4e);
        n4s = computeN4s(n4e);
        length4s = computeLength4s(c4n, n4s);
        
        % starting variables
        if level == 1
            uCR = zeros(nrSides, 1);
            Lambda4e = zeros(size(n4e,1), 1);
        else
            uCR = computeRefinementExtension(c4nOld,n4eOld,c4n,n4e,uCR);
            DuCR4e = computeGradientCR(uCR, gradsCR4e, s4e);
            Lambda4e = DuCR4e./repmat(sqrt(sum(DuCR4e.^2,2)),1,2);
            Lambda4e(isnan(Lambda4e)) = 0;
        end
        
        [uCR, Energy, nrSteps] = solvePrimalDualFormulation(uCR, Lambda4e, alpha, epsilon, tau, ...
                                STIMACR, MASSCR, intFCR14e, intFCR24e, intFCR34e,...
                                gradsCR4e, s4e, dof, nrSides, area4e);
        fprintf('Solution obtained for %i dof with minimal energy %.6f in %i iterations\n',...
                                                                nrDof, Energy, nrSteps);
        %% ESTIMATE ---------------------------------------------------------------------------
        [etaSq, etaSq4e, volError4e, jumpError4e] ...
                = computeEta(uCR, alpha, gamma, area4e, length4s, n4e, e4s, n4s,...
                             s4e, intFSq4e, intFCR14e, intFCR24e, intFCR34e);
        volError4lvl(level) = sum(volError4e);
        jumpError4lvl(level) = sum(jumpError4e);
        nrDof4lvl(level) = nrDof;
        eta4lvl(level) = sqrt(etaSq);
        energy4lvl(level) = Energy;
        
        % exact error
        if exactS
            errorL2Exact = error4eCRL2(c4n,n4e,uExact,uCR);
            errorL2Exact = sqrt(sum(errorL2Exact));
            errorL2Exact4lvl(level) = errorL2Exact;
            energyErr4lvl(level) = abs(Energy - EnergyExact);
        end
        
        % break condition
        if nrDof > minDof
            break;
        else
            level = level + 1;
        end
        
        %% MARK -------------------------------------------------------------------------------
        n4sMarked = markBulk(n4e, etaSq4e, theta);
        %% REFINE -----------------------------------------------------------------------------
        if theta == 1
            c4nOld = c4n;
            n4eOld = n4e;
            
            [c4n, n4e, n4sDb, ~] = refineUniformRed(c4n, n4e, n4sDb, zeros(0,2));
        else
            c4nOld = c4n;
            n4eOld = n4e;
            
            [c4n,n4e,n4sDb,~] = refineBi3GB(c4n,n4e,n4sDb,zeros(0,2),n4sMarked);
        end
    end
    
    %% SAVE
    if exist(folderName,'dir') == 0
        mkdir(folderName); 
    end
    
    % triangulation
    % c4n
    filename = [folderName,'/c4n.dat'];
    fid=fopen(filename,'wt');
    fprintf(fid,'%s\t%s\n','x','y');
    for j = 1:size(c4n,1)
        fprintf(fid,'%f\t%f\n',c4n(j,:));
    end
    fclose(fid);

    % n4e
    filename = [folderName,'/n4e.dat'];
    dlmwrite(filename,[n4e-1],'delimiter','\t','newline','unix');

    % s4e
    filename = [folderName,'/s4e.dat'];
    fid=fopen(filename,'wt');
    for j = 1:size(s4e,1)
        fprintf(fid,'%f\t%f\t%f\n',s4e(j,:));
    end
    fclose(fid);
    
    % solution
    if exactS
        Iu = uExact(c4n);
        filename = [folderName,'/uExact.dat'];
        dlmwrite(filename,[c4n(n4e',:),reshape(Iu(n4e'),3*size(n4e,1),1)],...
                 'delimiter','\t','newline','unix');
             
        filename = [folderName,'energyErr.dat'];
        fid=fopen(filename,'wt');
        fprintf(fid,'%s\t%s\n','x','y');
        for j = 1:level
            fprintf(fid,'%f\t%.12f\n',nrDof4lvl(j),energyErr4lvl(j));
        end
        fclose(fid);
    end

    % eta
    filename = [folderName,'/eta.dat'];
    fid=fopen(filename,'wt');
    fprintf(fid,'%s\t%s\n','x','y');
    for j = 1:level
        fprintf(fid,'%f\t%.12f\n',nrDof4lvl(j),eta4lvl(j));
    end
    fclose(fid);
    
    % energy
    filename = [folderName,'/energy.dat'];
    fid=fopen(filename,'wt');
    fprintf(fid,'%s\t%s\n','x','y');
    for j = 1:level
        fprintf(fid,'%f\t%.12f\n',nrDof4lvl(j),energy4lvl(j));
    end
    fclose(fid);
    
    %% PLOT
    % triangulation
    triangulation = figure;
    plotTriangulation (c4n, n4e);
    savefig(triangulation, [folderName,'\Triangulation.fig']);
    
    % solution uCR (component-wise)
    fig = figure;
    plotCR(c4n, n4e, uCR, 'discrete solution');
    name = '\discreteSol.fig';
    savefig(fig, [folderName, name]);
    
    % exact solution
    if exactS
        exact = figure;
        plotP1(c4n, n4e, Iu, "nodal interpolation of exact solution");
        savefig(exact,[folderName,'\exactSol.fig']);
    end
    
    % convergence plot
    convFigure = figure;
    plotConvergence(nrDof4lvl, eta4lvl, 'eta');
    hold on;
    plotConvergence(nrDof4lvl, volError4lvl, 'volume error');
    plotConvergence(nrDof4lvl, jumpError4lvl, 'jump error');
    if exactS
        plotConvergence(nrDof4lvl, errorL2Exact4lvl, 'L2 exact error');
        plotConvergence(nrDof4lvl, energyErr4lvl, 'energy error')
    end
    hold off;
    savefig(convFigure,[folderName,'\Convergence.fig']);
    
end