function signs4e=computeSigns4e(c4n,n4e,n4s,s4e)
%% computeSigns4e - normal signs for element
% This function returns a matrix in which each row contains for each
% triangle the 3 signs of the scalarproduct of unit outer normal and the
% unit normal vector of the side.
%
% See also computeNormal4s, computeNormal4s.
%

% Copyright 2014 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

normal4e=permute(computeNormal4e(c4n,n4e),[3 2 1]);
normal4s=computeNormal4s(c4n,n4s);

signs4e=[sum(normal4e(:,:,1).*normal4s(s4e(:,1),:),2),...
         sum(normal4e(:,:,2).*normal4s(s4e(:,2),:),2),...
         sum(normal4e(:,:,3).*normal4s(s4e(:,3),:),2)];

end
