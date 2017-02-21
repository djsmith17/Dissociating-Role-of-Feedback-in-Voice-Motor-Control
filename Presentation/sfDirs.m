function dirs = sfDirs(project, expType)
% function dirs = fsDirs(project)
%
% Setting the path to code, data in and out, and other toolboxes.
%
% INPUT:
% project:     string. Name of the project run.
%
% OUTPUT:
% dirs
% 
% Andres        :       v1      : init. 24 Jan 2017

if nargin == 0, project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; end

%% Determine hostname of system 
[~,host] = system('hostname');
host     = deblank(host);

%% Set appropriate directories for code, data input and output, based on system hostname.
switch host
    case 'tongue';
        dirs.Data           = fullfile('C:\Users\djsmith\Documents\Pilot Data', project, expType);         % Dir w/ raw datafiles
        dirs.Results        = fullfile('C:\Users\djsmith\Documents\Pilot Results', project, expType);   % Dir to output analyzed datafiles and figures to
        dirs.Prelim         = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Presentation\PrelimFiles\'); %Dir for project specific helper files
%         dirs.DataTestIn     = fullfile('C:\Users\salacho\Documents\Speech Lab',project,'data\test');    % Dir w/ raw datafiles
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB',project,'Analysis');           % Dir w/ data analysis Code
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        dirs.saveFileRoot   = '';                                                                           % Used to name figures. Gives the type of analysis done, channels, amongs others
        dirs.saveFileSuffix = '';                                                                           % Used to name figures. Includes the filter, window size, signal processing
    case ''
        % add your host name here
    otherwise
end

%% Set up path so code is accessible to Matlab
addpath(genpath(dirs.Code));        % Add dir w/ your code path
addpath(genpath(dirs.Prelim));      % Add dir w/ Prelim Files
addpath(genpath(dirs.helpers));     % Add dir w/ helpers code path

end