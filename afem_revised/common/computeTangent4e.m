function tangent4e = computeTangent4e(c4n,n4e)
%% computeTangent4e - Tangents for elements.
%   computeTangent4e(c4n, n4e) computes the tangent vectors of all sides of all
%                          elements in the decomposition. The tangent
%                          vectors are normed and point in counterclockwise
%                          direction. c4n and n4e are as specified in the
%                          documentation.
%
%   See also: computeTangent4s, computeNormal4e
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if isempty(n4e)
        tangent4e = zeros(0,2);
        return;
    end

    %% Compute tangent4e.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    c4start  = c4n(allSides(:,1),:);
    c4end    = c4n(allSides(:,2),:);
    lengths  = sqrt(sum((c4end-c4start).^2,2));
    tangents = (c4end - c4start)./[lengths lengths];
    tangent4e(1,:,:) = tangents(1:size(n4e,1),:)';
    tangent4e(2,:,:) = tangents(size(n4e,1)+1:2*size(n4e,1),:)';
    tangent4e(3,:,:) = tangents(2*size(n4e,1)+1:3*size(n4e,1),:)';
end