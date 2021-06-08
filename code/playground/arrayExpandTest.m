function arrayExpandTest(nrEntries, chunkSize)
  format long;
  disp(sprintf('\nExpand for every single entry:'))
  tic;
  a = [];
  for j = 1:nrEntries
    a(end+1) = j;
  end
  t1 = toc;
  disp(t1);

  disp(sprintf('\nPreallocate Memory chunkwise:'))
  tic;
  a = NaN([1,chunkSize]);
  for j = 1:nrEntries
    if j>length(a)
      a = [a,NaN([1,chunkSize])];
    end
  a(j) = j;
  end
  a = a(1:j);
  t2 = toc;
  disp(t2);

  disp('====================')

  fprintf('Time for chunkwise preallocating is %g percent of time for expanding\n', ...
    t2/t1*100);
end

