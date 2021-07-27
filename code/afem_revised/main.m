function S = main(B, kU, kQ, kV,  plotU, getBaseU, getBaseV, getBaseQ, identifier)
%MAIN executes AFEMRun_PoissonDPG with the given input or the default
%values if none given.
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
%                         solvePoissonDPG.m (default: @getP1ShapeFunctions)
%           getBaseV    - function handle which delivers baseV, gradV,
%                         nCPv, nIPv as described in solvePoissonDPG.m
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
%                         data on each level.he

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

%% INITIALIZATION
initAFEM;

% get default values if not enough input
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

if nargin < 1
    B = Poisson_SquareExact(1); 
end
%% EXECUTION
S = AFEMRun_PoissonDPG(B, kU, kQ, kV,plotU, getBaseU, getBaseV, getBaseQ, identifier);
end

