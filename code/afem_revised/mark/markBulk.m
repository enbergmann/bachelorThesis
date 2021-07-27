function n4sMarked = markBulk(n4p,eta4p,OPTtheta)
%% markBulk - Mark given parts using the bulk criterion.
%   n4sMarked = markBulk(n4p, eta4p, OPTtheta) marks the parts with the
%       largest estimated errors. The aggregate error of the marked parts
%       is OPTtheta times the overall error. By default, OPTtheta is set
%       to 0.5. n4p contains the nodes for the parts: 2 entries per row
%       for sides, 3 entries per row for elements.
%       The output is a list of marked sides given by their end nodes.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if (nargin < 3)
        theta = 0.5;
    else
        theta = OPTtheta;
    end
    dimParts = size(n4p,2);

    %% Bulk criterion
    [eta4p,ind] = sort(eta4p,'descend');
    % avoid round-off errors between sum and cumsum (esp. if theta=1)
    cumsumEta4p = cumsum(eta4p);
    J = find(cumsumEta4p >= theta*cumsumEta4p(end),1,'first');
    I = ind(1:J);

    %% Mark sides
    if dimParts == 2 % Sides were given. Mark 'em all.
        n4sMarked = n4p(I,:);
    elseif dimParts == 3 % Elements were given. Mark all sides of these.
        allSidesMarked = [n4p(I,[1 2]);n4p(I,[2 3]);n4p(I,[3 1])];
        % Eliminate duplicates.
        [b, ind]   = unique(sort(allSidesMarked,2), 'rows');
        n4sMarked = allSidesMarked(sort(ind),:);
    end
end
