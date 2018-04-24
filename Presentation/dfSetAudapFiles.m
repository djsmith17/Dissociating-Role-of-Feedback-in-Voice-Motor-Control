function audStimP = dfSetAudapFiles(expParam, dirs, trial, debug)
% audStimP = dfSetAudapFiles(expParam, dirs, trial, debug) sets the ost and
% pcf functions for a custom pitch-shift reflex experiment recorded in 
% Audapter. 
% This script creates two possible pitch-shifts. The first is a linear ramp
% stimulus that has a rise/fall time of 150ms, and magnitude of 100 cents.
% The second is sigmoidal ramp stimulus that matches the rise/fall time 
% and magnitude of the laryngeal perturbation experiment.

trialType = expParam.trialType(trial);
trigs     = expParam.trigs(trial,:,1);

ost       = expParam.ostFN;
pcf       = expParam.pcfFN;
trialLen  = expParam.trialLen;
pertName  = expParam.AudPert;
pertSw    = expParam.AudPertSw;

InflaT = expParam.InflaT;
InflaV = expParam.InflaV;

audStimP = organizeStimulus(trialType, trialLen, trigs, pertSw, pertName, InflaT, InflaV);

OST_tline = writeOSTportions(audStimP);
PCF_tline = writePCFportions(audStimP);

svPSRLevels(ost, OST_tline);
svPSRLevels(pcf, PCF_tline);

if debug
    drawStimulus(audStimP, dirs)
end
end

function audStimP = organizeStimulus(trialType, trialLen, trigs, pertSw, pertName, InflaT, InflaV)

audStimP.trialType = trialType; % 0: Control 1: Catch
audStimP.pertSw    = pertSw;    % 0: -100c   1: LarMag
audStimP.pertName  = pertName;  % 'Linear Standard', 'Sigmoid Matched'
audStimP.tStep     = 0.001;    % seconds
audStimP.fs        = 1/audStimP.tStep;
audStimP.lenTrialT = trialLen;                                  % Trial Length (Seconds)
audStimP.lenTrialP = audStimP.lenTrialT*audStimP.fs;            % Trial Length (Pert-Points)
audStimP.time      = (0:1:audStimP.lenTrialP-1)/audStimP.fs; % Seconds (Pert-Points)

audStimP.StTime    = trigs(1);                          % Seconds
audStimP.SpTime    = trigs(2);                          % Seconds
audStimP.StPoint   = round(audStimP.StTime*audStimP.fs);% Points
audStimP.SpPoint   = round(audStimP.SpTime*audStimP.fs);% Points

audStimP.InflaT    = InflaT;  % seconds
audStimP.InflaV    = InflaV;  % cents
audStimP.rampLenT  = [];      % Time (s)
audStimP.rampT     = [];      % Time Vector of ramp from Start (0) to finish (rampLenT)
audStimP.rampLenP  = [];      % Points
audStimP.rampMin   = [];
audStimP.ramp      = [];
audStimP.rampRv    = [];

audStimP.steadyLenP  = [];
audStimP.steadyLenT = [];

%Define the slope for the Aud. perturbation stimulus
if pertSw == 0 %Linear Standard Stimulus
    audStimP.rampLenT   = 0.15; % seconds          HardSet
    audStimP.rampT      = 0:audStimP.tStep:audStimP.rampLenT;
    audStimP.rampLenP   = length(audStimP.rampT);
    
    if trialType == 0
        audStimP.rampMin = 0;
        audStimP.ramp    = zeros(audStimP.rampLenP, 1);
    else
        audStimP.rampMin = -100; % cents            HardSet
        audStimP.ramp    = linspace(0, audStimP.rampMin, audStimP.rampLenP);
    end              
elseif pertSw == 1 %Sigmoid (Laryngeal) Matched Stimulus
    audStimP.rampLenT   = InflaT; % seconds
    audStimP.rampT      = 0:audStimP.tStep:audStimP.rampLenT;
    audStimP.rampLenP   = length(audStimP.rampT);
    
    if trialType == 0
        audStimP.rampMin = 0;
        audStimP.ramp    = zeros(audStimP.rampLenP, 1);
    else
        x = linspace(0, 10, audStimP.rampLenP);
        audStimP.rampMin  = InflaV; % cents
        audStimP.ramp     = audStimP.rampMin*sigmf(x, [1 5]);
    end                     
end

audStimP.ramp    = round((audStimP.ramp/100), 3);
audStimP.rampMin = round((audStimP.rampMin/100), 3);
audStimP.rampRv  = fliplr(audStimP.ramp);

audStimP.rampDNRange = audStimP.StPoint + (0:audStimP.rampLenP-1);
audStimP.rampUPRange = audStimP.SpPoint + (0:audStimP.rampLenP-1);
audStimP.steadySt    = audStimP.rampDNRange(end)+1;
audStimP.steadySp    = audStimP.rampUPRange(1)-1;
audStimP.steadyRange = audStimP.steadySt:audStimP.steadySp;
audStimP.steadyLenP  = length(audStimP.steadyRange);
audStimP.steadyLenT  = audStimP.steadyLenP/audStimP.fs;

stim = zeros(audStimP.lenTrialP,1);
stim(audStimP.rampDNRange) = audStimP.ramp;
stim(audStimP.steadyRange) = audStimP.rampMin;
stim(audStimP.rampUPRange) = audStimP.rampRv;

audStimP.stim = stim;
end

function OST_tline = writeOSTportions(audStimP)
%The Online Status Tracking file regulates the timing of when actions or
%changes to the speech occur. The steps between timed actions followed
%rules outlined in the Audapter Manuel. This has been specifically
%organized for customized Pitch-Shift Reflex experiments.

tStep      = audStimP.tStep;
StTime     = audStimP.StTime;
rampLen    = audStimP.rampLenP;
steadyLen  = audStimP.steadyLenT;

%The number of changes to f0 + the hold + last THREE clean-up lines
n = 2*rampLen + 1 + 3;

%p = pre-experiment lines in OST file
p = 7;

%The first 7 lines (p) of the OST should be nearly the same for all 
%participants. rmsSlopeWin might change eventually. The experiment starts 
%with a wait for voicing and then a 1.7-2.1s pause. 
OST_tline{1} = '# Online status tracking (OST) configuration file';
OST_tline{2} = 'rmsSlopeWin = 0.030000';
OST_tline{3} = ' ';
OST_tline{4} = '# Main section: Heuristic rules for tracking';
OST_tline{5} = ['n = ' num2str(n)];
OST_tline{6} = '0 INTENSITY_RISE_HOLD 0.01 0.05 {} # Detect voicing onset';
OST_tline{7} = ['2 ELAPSED_TIME ' num2str(StTime) ' NaN {} #The amount of time pre-perturbation']; %Random start between 1.7 and 2.1s

%The +2 comes from the numbering on the OST ahead of these commands
for i = 1:n
    if i <= rampLen
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(tStep) ' NaN {} #DownShift ' num2str(i) ' of ' num2str(rampLen)];
    elseif i == rampLen + 1 
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(steadyLen) ' NaN {} #Hold for the pitch-shift hold period'];
    elseif i <= 2*rampLen + 1
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(tStep) ' NaN {} #UpShift ' num2str(i) ' of ' num2str(rampLen)];    
    elseif i == 2*rampLen + 2
        OST_tline{i+p} = [num2str(i+2) ' OST_END NaN NaN {} #End the dang thing'];
    elseif 1 == 2*rampLen + 3
        OST_tline{i+p} = ' ';
    elseif i == 2*rampLen + 4
        OST_tline{i+p} = 'n = 0';
    end    
end
% Delete the below when total confidence in file creation
% OST_tline{8} = '3 ELAPSED_TIME 0.035 NaN {} #Shift 1 of 8';
% OST_tline{9} = '4 ELAPSED_TIME 0.035 NaN {} #Shift 2 of 8';
% OST_tline{10} = '5 ELAPSED_TIME 0.035 NaN {} #Shift 3 of 8';
% OST_tline{11} = '6 ELAPSED_TIME 0.035 NaN {} #Shift 4 of 8';
% OST_tline{12} = '7 ELAPSED_TIME 0.035 NaN {} #Shift 5 of 8';
% OST_tline{13} = '8 ELAPSED_TIME 0.035 NaN {} #Shift 6 of 8';
% OST_tline{14} = '9 ELAPSED_TIME 0.035 NaN {} #Shift 7 of 8';
% OST_tline{15} = '10 ELAPSED_TIME 0.035 NaN {} #Shift 8 of 8';
% OST_tline{16} = '11 ELAPSED_TIME 0.42 NaN {} #Hold for the trough';
% OST_tline{17} = '12 ELAPSED_TIME 0.035 NaN {} #Shift 1 of 8';
% OST_tline{18} = '13 ELAPSED_TIME 0.035 NaN {} #Shift 2 of 8';
% OST_tline{19} = '14 ELAPSED_TIME 0.035 NaN {} #Shift 3 of 8';
% OST_tline{20} = '15 ELAPSED_TIME 0.035 NaN {} #Shift 4 of 8';
% OST_tline{21} = '16 ELAPSED_TIME 0.035 NaN {} #Shift 5 of 8';
% OST_tline{22} = '17 ELAPSED_TIME 0.035 NaN {} #Shift 6 of 8';
% OST_tline{23} = '18 ELAPSED_TIME 0.035 NaN {} #Shift 7 of 8';
% OST_tline{24} = '19 ELAPSED_TIME 0.035 NaN {} #Shift 8 of 8';


% OST_tline{25} = '20 OST_END NaN NaN {} #End the dang thing';
% OST_tline{26} = ' ';
% OST_tline{27} = 'n = 0';
end

function PCF_tline = writePCFportions(audStimP)
%The Pertrubation Configuration file defines the levels for acoustic 
%variables at each action step defined in the OST. This have been
%specifically organized for customized Pitch-Shift Reflex experiments. 
%The PCF expects f0 in units of semitones. My analysis saves f0 in cents. 
%Divide by 100 to convert.

rampLen   = audStimP.rampLenP;
ramp      = audStimP.ramp;
rampMin   = audStimP.rampMin;
rampRv    = audStimP.rampRv;

%The number of changes to f0 + the hold + the last ONE clean-up line
n = 2*rampLen + 1 + 1;

%p = pre-experiment lines in PCF file
p = 8;

%The first 8 lines (p) of the PCF should be nearly the same for all 
%participants. No time warping is present and the experiment starts with a
%wait for voicing and then a 1.7-2.1s pause. 
PCF_tline{1} = '# Section 1 (Time warping): tBegin, rate1, dur1, durHold, rate2';
PCF_tline{2} = '0';
PCF_tline{3} = ' ';
PCF_tline{4} = '# Section 2: stat pitchShift(st) gainShift(dB) fmtPertAmp fmtPertPhi(rad)';
PCF_tline{5} = num2str(n+3);
PCF_tline{6} = '0, 0.0, 0.0, 0, 0';
PCF_tline{7} = '1, 0.0, 0.0, 0, 0';
PCF_tline{8} = '2, 0.0, 0.0, 0, 0';

%The +2 comes from the numbering on the PCF ahead of these commands
for i = 1:n
    if i <= rampLen
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(ramp(i)) ', 0.0, 0, 0'];
    elseif i == rampLen + 1 
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(rampMin) ', 0.0, 0, 0'];
    elseif i <= 2*rampLen + 1
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(rampRv(i-rampLen-1)) ', 0.0, 0, 0'];
    elseif i == 2*rampLen + 2
        PCF_tline{i+p} = [num2str(i+2) ', 0.0, 0.0, 0, 0'];
    end       
end
% Delete the below when total confidence in file creation
% PCF_tline{9} = '3, -0.03, 0.0, 0, 0';
% PCF_tline{10} = '4, -0.04, 0.0, 0, 0';
% PCF_tline{11} = '5, -0.04, 0.0, 0, 0';
% PCF_tline{12} = '6, -0.05, 0.0, 0, 0';
% PCF_tline{13} = '7, -0.18, 0.0, 0, 0';
% PCF_tline{14} = '8, -0.39, 0.0, 0, 0';
% PCF_tline{15} = '9, -0.55, 0.0, 0, 0';
% PCF_tline{16} = '10, -0.56, 0.0, 0, 0';
% PCF_tline{17} = '11, -0.56, 0.0, 0, 0';
% PCF_tline{18} = '12, -0.03, 0.0, 0, 0';
% PCF_tline{19} = '13, -0.04, 0.0, 0, 0';
% PCF_tline{20} = '14, -0.04, 0.0, 0, 0';
% PCF_tline{21} = '15, -0.05, 0.0, 0, 0';
% PCF_tline{22} = '16, -0.18, 0.0, 0, 0';
% PCF_tline{23} = '17, -0.39, 0.0, 0, 0';
% PCF_tline{24} = '18, -0.55, 0.0, 0, 0';
% PCF_tline{25} = '19, -0.56, 0.0, 0, 0';
% PCF_tline{26} = '20, 0.0, 0.0, 0, 0';
end

function tline = svPSRLevels (file, tline)

fid = fopen(file, 'w');
for i = 1:length(tline)    
    if i ~= length(tline)
        fprintf(fid,'%s\n',tline{i});
    else
        fprintf(fid,'%s',tline{i});
    end
end 
fclose(fid);

end

function drawStimulus(audStimP, dirs)
% close all

time     = audStimP.time;
stim     = audStimP.stim;

pertName = audStimP.pertName;
pertMag  = audStimP.rampMin;
pertTime = audStimP.rampLenT;

pertAx  = [audStimP.StTime, audStimP.SpTime];
pertAy  = [200 200];

plotpos = [10 50];
plotdim = [1200 500];
pertColor = [0.8 0.8 0.8];
AudStim = figure('Color', [1 1 1]);
set(AudStim, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on

plot(time, stim, 'LineWidth', 3)

xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('f0 Shift (cents)', 'FontSize', 16, 'FontWeight', 'bold')
title({'Auditory Feedback Perturbation Stimulus'; pertName}, 'FontSize', 18, 'FontWeight', 'bold')
axis([0 4 -1.1 1]); box off;

annotation('textbox',[0.62 0.25 0.40 0.1],...
           'String', {['Fall/Rise Time: ' num2str(pertTime) ' seconds'],...
                      ['Perturbation Magnitude: ' num2str(pertMag) ' cents']},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',14,...
                    'FontName','Arial');

set(gca, 'FontSize', 14,...
         'FontWeight', 'bold');
     
     
plTitle = ['AFPerturbStim_' pertName '.jpg'];     
saveFileName = fullfile(dirs.SavResultsDir, plTitle);
export_fig(saveFileName) 
end