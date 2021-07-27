function poolobj = initAFEM(nWorkers, identifier)
%%INITAFEM initializes the AFEM package
%   INITAFEM() initializes the AFEM package with serial computing.
%
%   poolobj = INITAFEM(nWorkers) enables parallel computing with nWorkers
%   cores and returns the pool object handle poolobj of parpool.
%
%   poolobj = INITAFEM(nWorkers, identifier) additionally sets an
%   identifier for the current session (default: current system time).

% Copyright 2016 Philipp Bringmann
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


%% PROCEED INPUT
if nargin < 1; nWorkers = 1; end
if nargin < 2; identifier = datestr(clock,'yymmdd_HHMMss'); end

%% INITIALIZATION
fprintf('\n*** INITIALIZE AFEM PACKAGE ***\n\n');

% load path
fprintf('  Add %s to path recursively ... ', pwd);
addpath(genpath(pwd));
fprintf('done\n');

% enable logging
fprintf('  Enable logging to logs/%s.log ... ', identifier);
diary(['logs/', identifier, '.log']);
fprintf('done\n');

% enable parallel computing
if nWorkers > 1
  fprintf('  Enable parallel computing ... \n');
  if ~exist('logs','dir'); mkdir('logs'); end
  if isempty(gcp('nocreate'))
    poolobj = parpool('local', nWorkers);
    fprintf('done\n');
  else
    fprintf('already active\n');
  end
elseif isempty(gcp('nocreate'))
  fprintf('  Use serial computing ... ok\n');
else
  fprintf('  Disable parallel computing ... ');
  %delete(gcp('nocreate'));
  fprintf('done\n');
end

fprintf('\n');

end
