function command = c

  if isunix()
    name = getenv('USER');
  else
    name = getenv('username');
  end

  command = '';
  if strcmp(name,'bergmaen') || strcmp(name,'Enrico') || strcmp(name,'enrico')
    if feature('IsDebugMode')
      % evalin('caller','dbquit all');
      command = 'dbquit all'
    end
    close all
    evalin('caller', 'clear all');
    clc
  end
end


