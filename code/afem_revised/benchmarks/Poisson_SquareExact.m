function B = Poisson_SquareExact(nRefinements)
% 
%
%


%% PROCEED INPUT
if nargin < 1; nRefinements = 1; end

%% INITIALIZATION
B = struct();

%% DEFAULT PARAMETERS
B.problem = 'Poisson';
B.name = 'SquareExact';
B.minNdof = 1e5;
B.theta = 1;

% Gauss-Newton parameters
B.maxNIter = 3;
B.tol1 = 0.001;
B.tol2 = 1;
B.rho = 0.5;
B.alphaMin = 0.001;

%% GEOMETRY
[B.c4n, B.n4e, B.n4sDb, B.n4sNb] = ...
  loadGeometry('Square', nRefinements);

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
%   val = 2*x(:,1) - 2*x(:,1).^2 + 2*x(:,2) - 2*x(:,2).^2;
     val = 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); %paper data
end
function val = u4Db(x)
  val = zeros(size(x, 1), 1);
%   val = sin(x(:,1))+sin(x(:,2));
end
function val = uExact(x)
%   val = x(:,1).*(1-x(:,2)).*x(:,2).*(1-x(:,1));
  val = sin(pi*x(:,1)).*sin(pi*x(:,2)); %paper data
end
function val = gradUExact(x)
%   val = [x(:,2) - x(:,2).^2 - 2*x(:,1).*x(:,2) + 2*x(:,1).*x(:,2).^2,...
%          x(:,1) - x(:,1).^2 - 2*x(:,1).*x(:,2) + 2*x(:,2).*x(:,1).^2];
    val = pi* [cos(pi*x(:,1)).*sin(pi*x(:,2)),...
        cos(pi*x(:,2)).*sin(pi*x(:,1))];
end
