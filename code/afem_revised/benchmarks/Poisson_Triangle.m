function B = Poisson_Triangle(nRefinements)
%
%
%


%% PROCEED INPUT
if nargin < 1; nRefinements = 0; end

%% INITIALIZATION
B = struct();

%% DEFAULT PARAMETERS
B.problem = 'Poisson';
B.name = 'Triangle';
B.minNdof = 1000;
B.theta = 0.5; 

% Gauss-Newton parameters
B.maxNIter = 3;
B.tol1 = 0.001;
B.tol2 = 1;
B.rho = 0.5;
B.alphaMin = 0.001;

%% GEOMETRY
%[B.c4n, B.n4e, B.n4sDb, B.n4sNb] = ...
%  loadGeometry('Square', nRefinements);
B.c4n = [1 0; 0 1; 0 0 ];
B.n4e = [1 2 3];
B.n4sDb = [2 3; 3 1; 1 2];
B.n4sNb = [];

%plotTriangulation(B.c4n,B.n4e);
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
%   val = 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); %paper data
end
function val = u4Db(x)
  val = zeros(size(x, 1), 1); %paper data
end
