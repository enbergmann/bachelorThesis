function plotRT0(c4n, n4e, u, p, OPTtitle4u, OPTtitle4p)
%% plotRT0 plots the solution u and its flux p
% INPUT:  c4n, n4e      - the triangulation
%         u             - the discrete solution
%         p             - the flux of the discrete solution
%         OPTtitle4u    - optional parameter, which defines the title
%                         for the solution u. The default value is empty.
%         OPTtitle4p    - optional parameter, which defines the title
%                         for the flux p. The default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.


mid4e = computeMid4e(c4n, n4e);

    if nargin >= 5
                 title4u = OPTtitle4u;
           else
                 title4u = '';
    end

figure;
plotP04e(c4n, n4e, u , title4u);

figure;
quiver2(mid4e(:,1),mid4e(:,2),p(:,1),p(:,2),'n=',0.1,'w=',[1 1]);

   if nargin == 6
               title(OPTtitle4p);
         else
               title('');
   end
end


