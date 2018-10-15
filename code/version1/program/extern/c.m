function c
% FULL CLEAR USE WITH CAUTION!!
% can only be used by user 'bethke'
% -clear all
% -close all
% -clc

name = getUserName ();

if strcmp(name,'bergmaen') || strcmp(name,'Enrico') || strcmp(name,'enrico')
evalin('caller','clear all');
evalin('caller','close all');
evalin('caller','clc');
if feature('IsDebugMode')
evalin('caller','dbquit all');
end
end

end


function name = getUserName ()
if isunix()
name = getenv('USER');
else
name = getenv('username');
end
end
