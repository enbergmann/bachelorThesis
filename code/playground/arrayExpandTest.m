function arrayExpandTest(nrEntries)
  disp(sprintf('\nExpand for every single entry:'))
  tic;

  a = [];
  for j = 1:nrEntries
    a(end+1) = j;
  end
  toc

  disp(sprintf('\nPreallocate Memory chunkwise:'))
  tic;
  a = NaN([1,1000]);
  for j = 1:nrEntries
    if j>length(a)
      a = [a,NaN([1,1000])];
    end
  a(j) = j;
  end
  a = a(1:j);
  toc

  disp('====================')
end

