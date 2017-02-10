function dirs = sfDirs(paradigm)
% function dirs = fsDirs(paradigm)
%
% Setting the path to code, data in and out, and other toolboxes.
%
% INPUT:
% paradigm:     string. Name of the paradigm run.
%
% OUTPUT:
% dirs
% 
% Andres        :       v1      : init. 24 Jan 2017

if nargin == 0, paradigm = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; end

%% Determine hostname of system 
[~,host] = system('hostname');
host     = deblank(host);

%% Set appropriate directories for code, data input and output, based on system hostname.
switch host
    case 'tongue';
        dirs.Data           = fullfile('C:\Users\djsmith\Documents\Pilot Data',paradigm,'Somatosensory Perturbation_Perceptual');         % Dir w/ raw datafiles
        dirs.Results        = fullfile('C:\Users\djsmith\Documents\Pilot Results',paradigm,'Somatosensory Perturbation_Perceptual');   % Dir to output analyzed datafiles and figures to
        dirs.Prelim         = fullfile('C:\Users\djsmith\Documents\MATLAB', paradigm, 'Presentation\PrelimFiles\');
        dirs.DataTestIn     = fullfile('C:\Users\salacho\Documents\Speech Lab',paradigm,'data\test');    % Dir w/ raw datafiles
        dirs.Code           = fullfile('C:\Users\djsmith\Documents\MATLAB',paradigm,'Analysis');           % Dir w/ data analysis Code
        dirs.chronux        = 'C:\Users\salacho\Documents\MATLAB\Chronux_matlab\chronux_code_2_00';         % Add path to chronux toolbox
        dirs.eeglab         = 'C:\Users\salacho\Documents\MATLAB\eeglab9_0_4_5s';                           % Folder to eeglab toolbox
        dirs.helpers        = 'C:\Users\djsmith\Documents\MATLAB\MATLAB-Toolboxes';                         % Dir to multiple function used for general analysis
        dirs.saveFileRoot   = '';                                                                           % Used to name figures. Gives the type of analysis done, channels, amongs others
        dirs.saveFileSuffix = '';                                                                           % Used to name figures. Includes the filter, window size, signal processing
        %dirs.baxter         = '';                                                                           % Dir w/ baxter data analysis code
        %dirs.asciiLoc       = 'C:\Users\salacho\Documents\Data\layout\gtec_montage.mat';
    case ''
        % add your host name here
    otherwise
end

%% Set up path so code is accessible to Matlab
addpath(genpath(dirs.Code));        % Add dir w/ your code path
% addpath(genpath(dirs.chronux));     % Add dir w/chronux toolbox
% addpath(genpath(dirs.eeglab));      % Add dir w/eeglab toolbox
addpath(genpath(dirs.helpers));     % Add dir w/ helpers code path

end