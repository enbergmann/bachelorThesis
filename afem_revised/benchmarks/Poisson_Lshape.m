function B = Poisson_Lshape(nRefinements)
%
%
%


%% PROCEED INPUT
if nargin < 1; nRefinements = 0; end

%% INITIALIZATION
B = struct();

%% DEFAULT PARAMETERS
B.problem = 'Poisson';
B.name = 'Lshape';
B.minNdof = 1e4;
B.theta = 0.5;

% Gauss-Newton parameters
B.maxNIter = 3;
B.tol1 = 0.001;
B.tol2 = 1;
B.rho = 0.5;
B.alphaMin = 0.001;

%% GEOMETRY
[B.c4n, B.n4e, B.n4sDb, B.n4sNb] = ...
  loadGeometry('LshapeCrisCros', nRefinements);

%% PROBLEM INPUT DATA
B.f = @f;
B.degreeF = 0;
B.u4Db = @u4Db;
B.degreeUDb = 0;

%% EXACT SOLUTION
B.exactKnown = false;

end


%% PROBLEM INPUT DATA
function val = f(x)
  val = ones(size(x, 1), 1);
end
function val = u4Db(x)
  val = zeros(size(x, 1), 1);
end
