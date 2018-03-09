function dirs = dfDirs(project)
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
    case 'SAR-D-635-1528'
        dirs.RecData        = fullfile('C:\DATA', project);  % Dir to save raw Data to
        dirs.SavData        = fullfile('W:\Experiments\', project);                  % Dir to open raw Data from
                                                                                     % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.Project        = fullfile('C:\GitHub', project);
        dirs.Prelim         = fullfile('C:\GitHub', project, 'Presentation\PrelimFiles'); % Dir for project specific helper files
        dirs.Code           = fullfile('C:\GitHub', project, 'Analysis');                 % Dir w/ data analysis Code
        dirs.Results        = fullfile('C:\GitHub', project, 'Results');                  % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\GitHub\MATLAB-Toolboxes';                               % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case 'dhcp-wifi-8021x-155-41-87-124.bu.edu'
        dirs.RecData        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/Pilot_Data');  % Dir to save raw Data to
        dirs.SavData        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/Pilot_Data');                  % Dir to open raw Data from
                                                                                     % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.Project        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study', project);
        dirs.Prelim         = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study', project, 'Presentation\PrelimFiles'); % Dir for project specific helper files
        dirs.Code           = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study', project, 'Analysis');                 % Dir w/ data analysis Code
        dirs.Results        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study', project, 'Results');                  % Dir to output analyzed datafiles and figures to
        dirs.helpers        = '/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/MATLAB-Toolboxes';                               % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case 'CNS-WS5'
        dirs.RecData        = fullfile('C:\Users\djsmith\My Documents\DATA', project); % Dir w/ raw datafiles
        dirs.SavData        = fullfile('W:\Experiments\', project);
       
        dirs.Project        = fullfile('C:\Users\djsmith\My Documents\MATLAB', project);
        dirs.Prelim         = fullfile('C:\Users\djsmith\My Documents\MATLAB', project, 'Presentation\PrelimFiles'); %Dir for project specific helper files
        dirs.Code           = fullfile('C:\Users\djsmith\My Documents\MATLAB', project, 'Analysis');           % Dir w/ data analysis Code
        dirs.Results        = fullfile('C:\Users\djsmith\My Documents\MATLAB', project, 'Results');            % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\My Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case 'tongue'
        dirs.RecData        = fullfile('C:\Users\djsmith\Documents\DATA', project); % Dir w/ raw datafiles
        dirs.SavData        = fullfile('W:\Experiments\', project);
       
        dirs.Project        = fullfile('C:\Users\djsmith\Documents\MATLAB', project);
        dirs.Prelim         = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Presentation\PrelimFiles'); %Dir for project specific helper files
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Analysis');           % Dir w/ data analysis Code
        dirs.Results        = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Results');            % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';                                                                           % Used to name figures. Includes the filter, window size, signal processing
    case 'DanteRig'
        dirs.RecData        = fullfile('E:\Documents\DATA', project); % Dir w/ raw datafiles  
        dirs.SavData        = fullfile('W:\Experiments\', project);
        
        dirs.Project        = fullfile('E:\Documents\MATLAB', project);
        dirs.Prelim         = fullfile('E:\Documents\MATLAB', project, 'Presentation\PrelimFiles'); %Dir for project specific helper files
        dirs.Code           = fullfile('E:\Documents\MATLAB', project, 'Analysis');           % Dir w/ data analysis Code
        dirs.Results        = fullfile('E:\Documents\MATLAB', project, 'Results');            % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'E:\Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';  
    case '677-GUE-WL-0001'
        dirs.RecData        = fullfile('C:\Users\djsmith\Documents\DATA', project); % Dir w/ raw datafiles  
        dirs.SavData        = fullfile('W:\Experiments\', project);
        
        dirs.Project        = fullfile('C:\Users\djsmith\Documents\MATLAB');
        dirs.Prelim         = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Presentation\PrelimFiles'); %Dir for project specific helper files
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Analysis');           % Dir w/ data analysis Code
        dirs.Results        = fullfile('C:\Users\djsmith\Documents\MATLAB', project, 'Results');            % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = ''; 
    otherwise
        fprintf('\nERROR: Please set directories for this host\n')
        return
end

%% Set up path so code is accessible to Matlab
addpath(genpath(dirs.Code));        % Add dir w/ your code path
addpath(genpath(dirs.Prelim));      % Add dir w/ Prelim Files
addpath(genpath(dirs.helpers));     % Add dir w/ helpers code path
end