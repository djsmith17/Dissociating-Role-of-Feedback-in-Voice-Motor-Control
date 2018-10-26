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
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('C:\DATA', project);        % Dir to save raw Data to
        dirs.SavData        = fullfile('W:\Experiments', project); % Dir to open raw Data from
      
        dirs.Code           = fullfile('C:\GitHub', project);             % The full code base
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'C:\GitHub\dfResults\Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\GitHub\MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
        
    case 'dhcp-wifi-8021x-155-41-27-125.bu.edu'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/Pilot_Data');  % Dir to save raw Data to
        dirs.SavData        = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/Pilot_Data');  % Dir to open raw Data from
        
        dirs.Code           = fullfile('/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study', project);    % The full code base
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = '/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/dfResults/Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = '/Users/elainekearney/Dropbox/1--PostDoc_2018-2020/PD_Study/MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case 'CNS-WS5'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('C:\Users\djsmith\Documents\DATA', project);   % Dir to save raw Data to
        dirs.SavData        = fullfile('W:\Experiments\', project);                   % Dir to open raw Data from
        
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB', project); % The full code base
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'C:\Users\djsmith\Documents\MATLAB\dfResults\Results';  % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';   % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case 'DanteRig'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('E:\Documents\DATA', project); % Dir to save raw Data to
        dirs.SavData        = fullfile('W:\Experiments\', project);   % Dir to open raw Data from
        
        dirs.Code           = fullfile('E:\Documents\MATLAB', project);
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'E:\Documents\MATLAB\dfResults\Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'E:\Documents\MATLAB\MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';  
    case '677-GUE-WL-0001'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('C:\Users\djsmith\Documents\DATA', project); % Dir to save raw Data to
        dirs.SavData        = fullfile('W:\Experiments', project);                  % Dir to open raw Data from
        
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB', project);
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'C:\Users\djsmith\Documents\MATLAB\dfResults\Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = '';
    case '677-gue-wl-0003'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('C:\DATA', project);        % Dir to save raw Data to  
        dirs.SavData        = fullfile('W:\Experiments', project); % Dir to open raw Data from
        
        dirs.Code           = fullfile('C:\GitHub', project);
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'C:\GitHub\dfResults\Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\GitHub\MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
        dirs.RecFileDir     = '';
        dirs.RecWaveDir     = '';
        dirs.SavFileDir     = '';
        dirs.SavResultsDir  = '';
       
        dirs.InflaRespFile  = '';
        dirs.saveFileSuffix = ''; 
        
    case '677-GUE-WD-0002'
        % RecData must be moved to SavData for backup and local disk space consolidation
        dirs.RecData        = fullfile('C:\DATA', project);        % Dir to save raw Data to  
        dirs.SavData        = fullfile('W:\Experiments', project); % Dir to open raw Data from
        
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\GitHub', project);
        dirs.Presentation   = fullfile(dirs.Code, 'Presentation');        % The scripts required for presentation
        dirs.Prelim         = fullfile(dirs.Presentation, 'PrelimFiles'); % Dir for presentation setup files
        dirs.Analysis       = fullfile(dirs.Code, 'Analysis');            % Dir w/ Code for data analysis
        
        dirs.Results        = 'C:\Users\djsmith\Documents\GitHub\dfResults\Results'; % Dir to output analyzed datafiles and figures to
        dirs.helpers        = 'C:\Users\djsmith\Documents\GitHub\MATLAB-Toolboxes';  % Dir to multiple function used for general analysis
        
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
addpath(genpath(dirs.Results));     % Where to find and save results
addpath(genpath(dirs.helpers));     % Add dir w/ helpers code path
end