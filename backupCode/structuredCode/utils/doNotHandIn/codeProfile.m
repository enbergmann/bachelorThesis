function codeProfile(strArr, expName, nrCalls)
  if nargin<2, expName = 'noNameGiven'; end
  if nargin<3, nrCalls = 500; end
  profile on
  for j = 1:nrCalls
    if mod(j, 10)==0, fprintf('%d ', j); end
    if mod(j,100)==0, fprintf('\n'); end
    for cmd = 1:length(strArr), evalin('caller', strArr(cmd)); end
  end
  profsave(profile('info'), sprintf('profiles/%s', expName))
  profile viewer
end

