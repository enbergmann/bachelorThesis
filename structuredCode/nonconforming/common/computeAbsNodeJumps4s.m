function absNodeJumps4s = computeAbsNodeJumps4s(n4e,e4s,nodeValues4e)
  absNodeJumps4s = zeros(size(e4s));
  for side = 1:size(e4s,1)
    tPlus = e4s(side,1);
    tMinus = e4s(side,2);
    if tMinus ~= 0
      pos = 1;
      for indPlus = 1:3
        indMinus = find(n4e(tMinus,:) == n4e(tPlus,indPlus));
        if not(isempty(indMinus))
          absNodeJumps4s(side,pos) = abs(nodeValues4e(tPlus,indPlus) - ...
                                         nodeValues4e(tMinus,indMinus));
          pos = pos + 1;
        end
      end
    end
  end
end
