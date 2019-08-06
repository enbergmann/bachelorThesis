function B = Poisson_LshapeExactUniform( nRefinements )


%% PROCEED INPUT
if nargin < 1; nRefinements = 0; end

%% INITIALIZATION
B = struct();

%% DEFAULT PARAMETERS
B.problem = 'Poisson';
B.name = 'LshapeExactUnfiform';
B.minNdof = 1e6;
B.theta = 1;

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

%% PROBLEM INPUT DATA
B.f = @f;
B.degreeF = 2;
B.u4Db = @u4Db;
B.degreeUDb = 0;

%% EXACT SOLUTION
B.exactKnown = true;
B.uExact = @uExact;
B.gradUExact = @gradUExact;

end


%% PROBLEM INPUT DATA
function val = f(x)
     val = -sqrt(x(:,1).^2 + x(:,2).^2).^(-1);
end
function val = u4Db(x)
  val = sqrt(x(:,1).^2 + x(:,2).^2); 
end
function val = uExact(x)
  val = sqrt(x(:,1).^2 + x(:,2).^2); 
end
function val = gradUExact(x)
    val =  [x(:,1)./sqrt(x(:,1).^2 + x(:,2).^2), x(:,2)./sqrt(x(:,1).^2 + x(:,2).^2)];
end