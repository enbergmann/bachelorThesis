function [osc,osc4e,mean4e] = oscillations(f,c4n,n4e,degree_f,nrrefinements,refinement)
%% oscillations
% calculates elementwise oscillations
% osc4e(T) = h_T ||f - fmean||_L^2(T)
%
% Input:     f            function handle for right-hand side f
%            c4n          coordinates for the nodes of the mesh
%            n4e          nodes for the elements of the mesh
%            n4sDb        the nodes of the sides in the Dirichlet boundary
%            n4sNb        the nodes of the sides in the Neumann boundary
%            degree_f     quadrature order for integrate
%           nrrefinements number of refinements before calculation
%            refinement   type of refinements
%                           default : @refineUniformRed
%                           also possible: @refineDual                         
%
% Output:    osc          total oscillations
%            osc4e        elementwise oscillations
%            mean4e       elementwise integral mean
%
% (C) 2009--2013 C. Merdon, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

if nargin < 4, degree_f = 4; end
if nargin < 5, nrrefinements = 0; end
if nargin < 6, refinement = @refineUniformRed; end

parents4e = 1:size(n4e,1);
for j=1:nrrefinements
    [c4n,n4e,~,~,parents4e_finer]=...
        refinement(c4n,n4e,zeros(0,2),zeros(0,2));
    parents4e = parents4e(parents4e_finer);
end 

area4e=computeArea4e(c4n,n4e);
n4s=computeN4s(n4e);
length4s=computeLength4s(c4n,n4s);
s4e=computeS4e(n4e);
if nargin(f) == 1
    f = @(n4p,Gpts4p,Gpts4ref)(f(Gpts4p));
end
mean4e=integrate(c4n,n4e,f,degree_f);
mean4e=mean4e./(area4e*ones(1,size(mean4e,2)));
osc4e=integrate(c4n, n4e,@(n4p,Gpts4p,Gpts4ref)(sum((f(n4p,Gpts4p,Gpts4ref)-mean4e).^2,2)),2*degree_f);
diam4e=max(length4s(s4e),[],2);
osc4e=osc4e.*diam4e.^2;
osc4e=accumarray(parents4e(:),osc4e(:));
osc=sqrt(sum(osc4e));
end
