function B = Poisson_SlitExact(nRefinements)
%
%
%


%% PROCEED INPUT
if nargin < 1; nRefinements = 1; end

%% INITIALIZATION
B = struct();

%% DEFAULT PARAMETERS
B.problem = 'Poisson';
B.name = 'SlitExact';
B.minNdof = 1e5;
B.theta = 0.5;

% Gauss-Newton parameters
B.maxNIter = 3;
B.tol1 = 0.001;
B.tol2 = 1;
B.rho = 0.5;
B.alphaMin = 0.001;

%% GEOMETRY
[B.c4n, B.n4e, B.n4sDb, B.n4sNb] = ...
  loadGeometry('Slit', nRefinements);

%% PROBLEM INPUT DATA
B.f = @f;
B.degreeF = 0;
B.u4Db = @u4Db;
B.degreeUDb = 10;

%% EXACT SOLUTION
B.exactKnown = true;
B.uExact = @uExact;
B.gradUExact = @gradUExact;

end


%% PROBLEM INPUT DATA
function val = f(x)
  val = zeros(size(x,1),1);
end
function val = u4Db(x)
  val = uExact(x);
end
function val = uExact(x)
  [phi, r] = cart2pol(x(:,1), x(:,2));
  phi(phi<-eps) = phi(phi<-eps) + 2*pi;
  val = r.^(1/4).*sin(1/4*phi);
end
function val = gradUExact(x)
  [phi, r] = cart2pol(x(:,1), x(:,2));
  phi(phi<-eps) = phi(phi<-eps) + 2*pi;
  if ~isempty(phi(phi<-eps)) | ~isempty(phi(phi>2*pi))
    error(['Evaluation of gradUExact: ',...
           'Transformation to polar coordinates failed'])
  end
  val = zeros(size(x, 1), 2);
  val(:,1) = 1/4*r.^(-3/4).*cos(1/4*phi);
  val(:,2) = 1/4*r.^(-3/4).*sin(1/4*phi);
  val = [sum(val.*[-sin(phi), cos(phi)], 2),...
         sum(val.*[ cos(phi), sin(phi)], 2)];
  val(r==0,2) = Inf;
end
