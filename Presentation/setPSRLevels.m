function audStimP = setPSRLevels(route, tStep, ost, pcf, trialType, spans)
%This function will take care of the ost and the pcf function for a custom
%pitch-shift reflex experiment based off a previously recorded 'route' the
%participant's pitch takes when they're larynx is physically perturbed
%Last Updated: 12/21/2016
%12/21/2016: Added a conditional to check if it is a catch trial or
%control. Extended pre-perturbation period from 0.5-->2.0s
%01/19/2017: The onset of perturbation is not semi-random to start between
%1.7 and 2.1 s following phonation

audStimP = organizeStimulus(route, tStep, trialType, spans);

OST_tline = writeOSTportions(audStimP);
PCF_tline = writePCFportions(audStimP);

svPSRLevels(ost, OST_tline);
svPSRLevels(pcf, PCF_tline);

% drawStimulus(audStimP)
end

function audStimP = organizeStimulus(route, tStep, trialType, spans)

routeStps = length(route);
if trialType == 0;
    route = zeros(1, numpts);
end

audStimP.AudFs     = 48000;  %Hardset
audStimP.lenTrialT = 4;      %Trial Length (Seconds) %Hardset
audStimP.route     = route;
audStimP.tStep     = tStep;                %Length of time-step (Seconds)
audStimP.tStepP    = round(tStep*audStimP.AudFs); %Length of time-step (Points)
audStimP.routeStps = routeStps;               %How many time-steps
audStimP.trialType = trialType;            %0 for control %1 for Catch
audStimP.lenTrialP = audStimP.lenTrialT*audStimP.AudFs; %Trial Length (Points)
audStimP.StTime    = spans(1);                               %Seconds
audStimP.SpTime    = spans(2);                               %Seconds
audStimP.StPoint   = round(audStimP.StTime*audStimP.AudFs);  %Points
audStimP.SpPoint   = round(audStimP.SpTime*audStimP.AudFs);  %Points
audStimP.lenPerT   = audStimP.SpTime - audStimP.StTime;      %Seconds
audStimP.lenPerP   = round(audStimP.lenPerT*audStimP.AudFs); %Points

audStimP.time     = (0:1:audStimP.lenTrialP-1)/audStimP.AudFs; %Projected recorded time course (Points)

audStimP.routeStpsT   = audStimP.routeStps*audStimP.tStep; %How long the route takes (Seconds)
audStimP.routeStpsP   = audStimP.routeStpsT*audStimP.AudFs; %How long the route takes (Points)
audStimP.routeInSp  = audStimP.StPoint + audStimP.routeStpsP; %Point when the shift 'bottoms out'
audStimP.routeOutSt = audStimP.SpPoint - audStimP.routeStpsP; %Point when the shift begins to increase to baseline
audStimP.lenPerVall = audStimP.routeOutSt - audStimP.routeInSp; %Points between either pertrubation 'route' (Valley)

stim = zeros(1, audStimP.lenTrialP);

for i = 1:audStimP.routeStps
    for j = 1:(audStimP.tStepP)
        stim(audStimP.StPoint + j + (i-1)*audStimP.tStepP) = audStimP.route(i);
        stim(audStimP.routeOutSt + j +(i-1)*audStimP.tStepP) = audStimP.route(routeStps - (i-1));        
    end
end
stim(audStimP.routeInSp:audStimP.routeOutSt) = audStimP.route(routeStps);

audStimP.stim = stim;
end

function OST_tline = writeOSTportions(audStimP)
%The Online Status Tracking file regulates the timing of when actions or
%changes to the speech occur. The steps between timed actions followed
%rules outlined in the Audapter Manuel. This has been specifically
%organized for customized Pitch-Shift Reflex experiments. Due to the nature
%of how the route is calculated, the elapsed time for each action is 
%currently fixed.

finaTim = num2str(0.42);  %HardSet **Fix this next**

%The number of changes to f0 + the last FOUR clean-up lines
n = audStimP.routeStps + 4;

%p = pre-experiment lines in OST file
p = 7;

%The first 7 lines (p) of the OST should be nearly the same for all 
%participants. rmsSlopeWin might change eventually. The experiment starts 
%with a wait for voicing and then a 0.5s pause. 
OST_tline{1} = '# Online status tracking (OST) configuration file';
OST_tline{2} = 'rmsSlopeWin = 0.030000';
OST_tline{3} = ' ';
OST_tline{4} = '# Main section: Heuristic rules for tracking';
OST_tline{5} = ['n = ' num2str(n)];
OST_tline{6} = '0 INTENSITY_RISE_HOLD 0.01 0.05 {} # Detect voicing onset';
OST_tline{7} = ['2 ELAPSED_TIME ' num2str(audStimP.StTime) ' NaN {} #The amount of time pre-perturbation']; %Random start between 1.7 and 2.1s

for i = 1:n
    if i <= audStimP.routeStps
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(audStimP.tStep) ' NaN {} #Shift ' num2str(i) ' of ' num2str(audStimP.routeStps)];
    elseif i == audStimP.routeStps + 1 
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' finaTim ' NaN {} #Hold for the rest of the perturbation'];
    elseif i == audStimP.routeStps + 2
        OST_tline{i+p} = [num2str(i+2) ' OST_END NaN NaN {} #End the dang thing'];
    elseif 1 == audStimP.routeStps + 3
        OST_tline{i+p} = ' ';
    elseif i == audStimP.routeStps + 4
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
% OST_tline{16} = '11 ELAPSED_TIME 0.42 NaN {} #Hold for the rest of the perturbation';
% OST_tline{17} = '12 OST_END NaN NaN {} #End the dang thing';
% OST_tline{18} = ' ';
% OST_tline{19} = 'n = 0';
end

function PCF_tline = writePCFportions(audStimP)
%The Pertrubation Configuration file defines the levels for acoustic 
%variables at each action step defined in the OST. This have been
%specifically organized for customized Pitch-Shift Reflex experiments. 
%The change in f0 is in semitones. The route saves f0 in cents. Divide by
%100 to convert

%The number of changes to f0 + the last TWO clean-up lines
n = audStimP.routeStps + 2;

%p = pre-experiment lines in PCF file
p = 8;

%The first 8 lines (p) of the PCF should be nearly the same for all 
%participants. No time warping is present and the experiment starts with a
%wait for voicing and then a 0.5s pause. 
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
    if i <= audStimP.routeStps
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(audStimP.route(i)/100) ', 0.0, 0, 0'];
    elseif i == audStimP.routeStps + 1 
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(audStimP.route(end)/100) ', 0.0, 0, 0'];
    elseif i == audStimP.routeStps + 2
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
% PCF_tline{18} = '12, 0.0, 0.0, 0, 0';
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

function drawStimulus(audStimP)
close all
plotpos = [200 100];
plotdim = [1300 500];
AudStim = figure('Color', [1 1 1]);
set(AudStim, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

plot(audStimP.time, audStimP.stim)
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Fundamental Frequency Shift (st)', 'FontSize', 12, 'FontWeight', 'bold')
title('Pitch-Shift Reflex Experiment Stimulus', 'FontSize', 16, 'FontWeight', 'bold')
axis([0 4 -60 0]); box off;

set(gca, 'FontSize', 16,...
         'FontWeight', 'bold');

end