function S = AFEMRun_PoissonDPG(B, kU, kQ, kV,  plotU, getBaseU, getBaseV, getBaseQ, identifier)
%%AFEMRUN_POISSONDPG Script for Adaptive LS-FEM for Stokes equations
%   S = AFEMRUN_POISSONDPG(B, kU, kQ, kV, getBaseU, getBaseV, plotU, identifier)
%   runs the AFEM loop and returns a structure handle S with input data and
%   the output data on each level.
%
% input:    B           - input structure for the problem
%           kU          - total polynomial degree of P_{k_u}(\mathcal{T})
%                         (default: 1)
%           kQ          - polynomial degree for P_{k_q}(\mathcal{E}) as
%                         a part of the trial space (default: 0)
%           kV          - total polynomial degree of P_{k_v}(\mathcal{T})
%                         (default: 1)
%           getBaseU    - function handle which delivers baseU, gradU,
%                         nCPu, nIPu, c4BaseU as described in
%                         solveEstPoissonDPG.m (default: @getP1ShapeFunctions)
%           getBaseV    - function handle which delivers baseV, gradV,
%                         nCPv, nIPv as described in solveEstPoissonDPG.m
%                         (default: @getP1ShapeFunctions)
%           getBaseQ    - function handle which delivers baseQ as described
%                         in solveEstPoissonDPG.m (default:
%                         @getP0TraceFunctions)
%           plotU       - function handle to plot the solution of the
%                         primal dPG method for the Poisson model problem
%           identifier  - identifier string for the current computation
%                         (default: current system time)
%
% output:   S           - structure handle with input data and the output
%                         data on each level.

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

%% PROCEED INPUT
if nargin < 9
    identifier = datestr(now,30);
end

if nargin < 2
    kU = 1;
    kQ = 0;
    kV = 1;
end

if nargin < 8
    if kQ == 0
        getBaseQ = @getP0TraceFunctions;
    elseif kQ == 1
        getBaseQ = @getP1TraceFunctions;
    elseif kQ == 2
        getBaseQ = @getP2TraceFunctions;
    else
        fprintf('No basis provided.');
    end
end


if nargin < 6
    if kU == 1
        getBaseU = @getP1ShapeFunctions;
    elseif kU == 2
        getBaseU = @getP2ShapeFunctions;
    elseif kU == 3
        getBaseU = @getP3ShapeFunctions;
    else
        fprintf('No basis provided.');
    end
    if kV == 1
        getBaseV = @getP1ShapeFunctions;
    elseif kV == 2
        getBaseV = @getP2ShapeFunctions;
    elseif kV == 3
        getBaseV = @getP3ShapeFunctions;
    else
        fprintf('No basis provided.');
    end
end

if nargin < 5
    if kU == 1
        plotU = @plotP1;
    elseif kU == 2
        plotU = @plotP2;
    elseif kU == 3
        plotU = @plotP3;
    else
        fprintf('No plot function provided.');
    end
end
[baseU, gradU, nCPu, nIPu, c4BaseU] = getBaseU();
[baseV, gradV, nCPv, nIPv]       = getBaseV();
[baseQ] = getBaseQ();

%% INITIAL NOTICE
fprintf(['*** ADAPTIVE COURANT FEM FOR POISSON MODEL PROBLEM ***\n', ...
    'This is AFEMRun_PoissonDPG function with parameters\n\n', ...
    '  name       = %s\n\n', ...
    '  minNdof    = %d\n  theta      = %g\n\n', ...
    '  identifier = %s\n\n'], ...
    B.name, B.minNdof, B.theta, identifier);

%% INITIALIZATION
% initialize data structure
S = B;
S.identifier = identifier;
lvl = 0;

% triangulation
c4n = B.c4n;
n4e = B.n4e;
n4sDb = B.n4sDb;

% measure time
t = clock;

ndof4lvl = [];
eta4lvl = [];
errL24lvl = [];
errH14lvl = [];
errEnergy4lvl = [];

%% AFEM LOOP
while(true)
    %% UPDATE LEVEL
    fprintf('\nLEVEL %d (%s)', lvl, datestr(now));
    lvl = lvl + 1;
    
    %% SOLVE
    fprintf('\n  SOLVE');
    % solve and estimate COURANT FEM
    tic;
    [u, ndof, u4e, eta, eta4e] = solveEstPoissonDPG(B.f, B.u4Db, c4n, n4e, n4sDb, B.degreeF, ...
        baseU, gradU, nCPu, nIPu, c4BaseU, baseQ, baseV, gradV, nCPv, nIPv, kU, kQ, kV);
    runtime = toc;
    ndof4lvl(lvl) = ndof;
    fprintf('\n    ndof  = %d\n    time  = %.4f sec.', ndof, runtime); %break;
    
    %% ESTIMATE
    fprintf('\n  ESTIMATE');
    eta4lvl(lvl) = eta;
    fprintf('\n    eta   = %0.8g', eta);
    
    %% ERROR COMPUTATION
    if B.exactKnown
        fprintf('\n  ERROR')
        errL2 = errorL2dPG(B.uExact, u4e, c4n, n4e, kU, baseU, B.degreeF);
        errH1 = errorH1dPG(B.uExact, B.gradUExact, u4e, c4n, n4e, kU, baseU, gradU, B.degreeF);
        errEnergy = errorEnergyDPG(B.gradUExact, u4e, c4n, n4e, kU, baseU, gradU, B.degreeF);
        fprintf('\n    errL2   = %.8g',errL2);
        fprintf('\n    errH1   = %.8g',errH1);        
        fprintf('\n    errEnergy   = %.8g',errEnergy);
        if lvl ~= 1
            rateL2 = log2(errL24lvl(lvl-1)/errL2);
            rateH1 = log2(errH14lvl(lvl-1)/errH1);
            fprintf('\n    rateL2  = %.8g',rateL2);
            fprintf('\n    rateH1  = %.8g',rateH1);
        else
            rateL2 = NaN;
            rateH1 = NaN;
        end
    else
        errL2 = NaN;
        rateL2 = NaN;
        errH1 = NaN;
        rateH1 = NaN;
        errEnergy = NaN;
    end
    errL24lvl(lvl) = errL2;
    errH14lvl(lvl) = errH1;
    errEnergy4lvl(lvl) = errEnergy;
    %% SAVE DATA
    S.level(lvl) =...
        struct('c4n',    c4n,  'n4e',     n4e,     'n4sDb',  n4sDb, ...
        'ndof',   ndof, 'runtime', runtime, ...
        'u',       u,   'eta',    eta, 'errH1', errH1, 'rateH1', rateH1, ...
        'errL2',    errL2, 'rateL2', rateL2, 'errEnergy', errEnergy);
    % save data to file
    % AFEMSave(S);
    
    %% BREAK CONDITION
    if ndof >= B.minNdof
        fprintf('\n\nBREAK AS NDOF >= MINNDOF\n\n'); break;
    end
    
    %% MARK AND REFINE
    if B.theta == 1
        fprintf('\n  RED-REFINE\n');
        [c4n, n4e, n4sDb] = refineUniformRed(c4n, n4e, n4sDb, zeros(0,2));
    else
        fprintf('\n  MARK');
        n4sMarked = markBulk(n4e, eta4e, B.theta);
        fprintf('\n  CLOSURE');
        n4sMarked = closure(n4e, n4sMarked);
        fprintf('\n  Bi3GB-REFINE\n');
        [c4n, n4e, n4sDb] = ...
            refineBi3GB(c4n, n4e, n4sDb, zeros(0, 2), n4sMarked);
    end
end

fprintf('\nTotal time evolved : %.4f sec.\n',etime(clock,t));

%% POST-PROCESSING
fprintf('\n');

%% PLOTS
if isa(plotU, 'function_handle')
    if strcmp(input(':: Plot triangulation? [y/N] ', 's'), 'y')
        figure; plotTriangulation(c4n, n4e);
    end
    if strcmp(input(':: Plot solution? [y/N] ', 's'), 'y')
        n4s = computeN4s(n4e);
        figure; plotU(c4n, n4e, u(1:end-(kQ+1)*size(n4s,1)));
    end
    if strcmp(input(':: Plot convergence history? [y/N] ', 's'), 'y')
        figure; plotConvergence(ndof4lvl, eta4lvl, 'eta');
        if B.exactKnown
            figure; plotConvergence(ndof4lvl, errL24lvl, 'L2 error');
            figure; plotConvergence(ndof4lvl, errH14lvl, 'H1 error');
            figure; plotConvergence(ndof4lvl, errEnergy4lvl, 'Energy error');
        end
    end
end
end

