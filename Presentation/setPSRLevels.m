function audStimP = setPSRLevels(InflaRespRoute, tStep, ost, pcf, trialType, trigs, stimType)
%This function will take care of the ost and the pcf function for a custom
%pitch-shift reflex experiment recorded in Audapter. The custom 
%perturbation shape and magnitude is based off a previously recorded 
%'route' the participant's pitch takes when they're larynx is physically 
%perturbed.

audStimP = organizeStimulus(trialType, stimType, trigs, InflaRespRoute, tStep);

OST_tline = writeOSTportions(audStimP);
PCF_tline = writePCFportions(audStimP);

svPSRLevels(ost, OST_tline);
svPSRLevels(pcf, PCF_tline);

drawStimulus(audStimP)
end

function audStimP = organizeStimulus(trialType, stimType, trigs, InflaRespRoute, tStep)

audStimP.trialType = trialType; % 0: Control 1: Catch
audStimP.AudFs     = 48000;     % Hardset
audStimP.lenTrialT = 4;                                         %Trial Length (Seconds) %Hardset
audStimP.lenTrialP = audStimP.lenTrialT*audStimP.AudFs;         %Trial Length (Points)
audStimP.time      = (0:1:audStimP.lenTrialP-1)/audStimP.AudFs; %Projected recorded time course (Points)

audStimP.StTime    = trigs(1);                               %Seconds
audStimP.SpTime    = trigs(2);                               %Seconds
audStimP.StPoint   = round(audStimP.StTime*audStimP.AudFs);  %Points
audStimP.SpPoint   = round(audStimP.SpTime*audStimP.AudFs);  %Points
audStimP.lenPerT   = audStimP.SpTime - audStimP.StTime;      %Seconds
audStimP.lenPerP   = round(audStimP.lenPerT*audStimP.AudFs); %Points

audStimP.route     = InflaRespRoute;
audStimP.routetStep= tStep;
audStimP.routeStps = length(InflaRespRoute);
audStimP.routeLenT = audStimP.routetStep*audStimP.routeStps;

if stimType == 1
    if trialType == 0;
        audStimP.ramp = zeros(size(InflaRespRoute));
    else
        audStimP.ramp = InflaRespRoute;
    end
    
    audStimP.rampFin  = audStimP.ramp(end);
    audStimP.rampStps = length(audStimP.ramp);
    audStimP.tStep    = tStep;                       %Length of time-step (Seconds)
    audStimP.tStepP   = round(audStimP.tStep*audStimP.AudFs); %Length of time-step (Points)
    
    audStimP.rampLenT = audStimP.rampStps*audStimP.tStep;       %How long is the ramp (Seconds)
    audStimP.rampLenP = round(audStimP.rampLenT*audStimP.AudFs); %How long is the ramp (Points) 
    
    audStimP.rampDNsp = round(audStimP.StPoint + audStimP.rampLenP); %Point when ramp down ends (starts to bottom out)
    audStimP.rampUPst = round(audStimP.SpPoint - audStimP.rampLenP); %Point when ramp up starts (begins to climb towards baseline)
    
    audStimP.rampDNRange = []; audStimP.rampUPRange = [];
    audStimP.rampDNValues = []; audStimP.rampUPValues = [];
    
    for r = 1:audStimP.rampStps
        startDN = audStimP.StPoint + (r-1)*audStimP.tStepP;
        audStimP.rampDNRange(:,r) = startDN : startDN + audStimP.tStepP;
        audStimP.rampDNValues(r)  = audStimP.ramp(r);
        
        startUP = audStimP.rampUPst + (r-1)*audStimP.tStepP;
        audStimP.rampUPRange(:,r) = startUP : startUP + audStimP.tStepP;
        audStimP.rampUPValues(r)  = audStimP.ramp(audStimP.rampStps - (r-1));
    end       
    
    audStimP.lenPerVallP = audStimP.rampUPst - audStimP.rampDNsp; %Points between either pertrubation 'route' (Valley)
    audStimP.lenPerVallT = audStimP.lenPerVallP/audStimP.AudFs;   %Seconds between either pertrubation 'route' (Valley)
    audStimP.VallRange   = audStimP.rampDNsp:audStimP.rampUPst;
    audStimP.VallValue   = audStimP.rampFin;    
    
elseif stimType == 2 %0 - minPoint Sigmoid    
    if trialType == 0;
        audStimP.ramp = zeros(100, 1);
    else
        x = 0:0.1:10;
        audStimP.ramp = InflaRespRoute(end)*sigmf(x,[1 5]);
    end
        
    audStimP.rampFin  = audStimP.ramp(end);
    audStimP.rampStps = length(audStimP.ramp);
    audStimP.tStep    = audStimP.routeLenT/100;
    audStimP.tStepP   = round(audStimP.tStep*audStimP.AudFs);        
   
    audStimP.rampLenT = audStimP.rampStps*audStimP.tStep; %How long is the ramp (Seconds)
    audStimP.rampLenP = round(audStimP.rampLenT*audStimP.AudFs); %How long is the ramp (Points) 

    audStimP.rampDNsp = round(audStimP.StPoint + audStimP.rampLenP); %Point when ramp down ends (starts to bottom out)
    audStimP.rampUPst = round(audStimP.SpPoint - audStimP.rampLenP); %Point when ramp up starts (begins to climb towards baseline)
    
    audStimP.rampDNRange = []; audStimP.rampUPRange = [];
    audStimP.rampDNValues = []; audStimP.rampUPValues = [];
    
    for r = 1:audStimP.rampStps
        startDN = audStimP.StPoint + (r-1)*audStimP.tStepP;
        audStimP.rampDNRange(:,r) = startDN : startDN + audStimP.tStepP;
        audStimP.rampDNValues(r)  = audStimP.ramp(r);
        
        startUP = audStimP.rampUPst + (r-1)*audStimP.tStepP;
        audStimP.rampUPRange(:,r) = startUP : startUP + audStimP.tStepP;
        audStimP.rampUPValues(r)  = audStimP.ramp(audStimP.rampStps - (r-1));
    end 
    
    audStimP.lenPerVallP = audStimP.rampUPst - audStimP.rampDNsp; %Points between either pertrubation 'route' (Valley)
    audStimP.lenPerVallT = audStimP.lenPerVallP/audStimP.AudFs;   %Seconds between either pertrubation 'route' (Valley)
    audStimP.VallRange   = audStimP.rampDNsp:audStimP.rampUPst;
    audStimP.VallValue   = audStimP.rampFin; 
    
elseif stimType == 3 %0 - 100 linearly   
    if trialType == 0;
        audStimP.ramp = zeros(99, 1);
    else
        audStimP.ramp = -1*(0:100-1);
    end
        
    audStimP.rampFin  = audStimP.ramp(end);
    audStimP.rampStps = length(audStimP.ramp);
    audStimP.tStep    = audStimP.routeLenT/100;
    audStimP.tStepP   = round(audStimP.tStep*audStimP.AudFs);        
    
    audStimP.rampLenT = audStimP.rampStps*audStimP.tStep; %How long is the ramp (Seconds)
    audStimP.rampLenP = round(audStimP.rampLenT*audStimP.AudFs); %How long is the ramp (Points) 

    audStimP.rampDNsp = round(audStimP.StPoint + audStimP.rampLenP); %Point when ramp down ends (starts to bottom out)
    audStimP.rampUPst = round(audStimP.SpPoint - audStimP.rampLenP); %Point when ramp up starts (begins to climb towards baseline)
    
    audStimP.rampDNRange = []; audStimP.rampUPRange = [];
    audStimP.rampDNValues = []; audStimP.rampUPValues = [];

    for r = 1:audStimP.rampStps
        startDN = audStimP.StPoint + (r-1)*audStimP.tStepP;
        audStimP.rampDNRange(:,r) = startDN : startDN + audStimP.tStepP;
        audStimP.rampDNValues(r)  = audStimP.ramp(r);
        
        startUP = audStimP.rampUPst + (r-1)*audStimP.tStepP;
        audStimP.rampUPRange(:,r) = startUP : startUP + audStimP.tStepP;
        audStimP.rampUPValues(r)  = audStimP.ramp(audStimP.rampStps - (r-1));
    end 
    
    audStimP.lenPerVallP = audStimP.rampUPst - audStimP.rampDNsp; %Points between either pertrubation 'route' (Valley)
    audStimP.lenPerVallT = audStimP.lenPerVallP/audStimP.AudFs;   %Seconds between either pertrubation 'route' (Valley)
    audStimP.VallRange   = audStimP.rampDNsp:audStimP.rampUPst;
    audStimP.VallValue   = audStimP.rampFin; 
end

stim = zeros(1, audStimP.lenTrialP);
for i = 1:audStimP.rampStps
    stim(audStimP.rampDNRange(:,i)) = audStimP.rampDNValues(i);
    stim(audStimP.rampUPRange(:,i)) = audStimP.rampUPValues(i);        
end
stim(audStimP.VallRange) = audStimP.VallValue;

audStimP.stim = stim;
end

function OST_tline = writeOSTportions(audStimP)
%The Online Status Tracking file regulates the timing of when actions or
%changes to the speech occur. The steps between timed actions followed
%rules outlined in the Audapter Manuel. This has been specifically
%organized for customized Pitch-Shift Reflex experiments.

%The number of changes to f0 + the hold + last THREE clean-up lines
n = 2*audStimP.rampStps + 1 + 3;

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
OST_tline{7} = ['2 ELAPSED_TIME ' num2str(audStimP.StTime) ' NaN {} #The amount of time pre-perturbation']; %Random start between 1.7 and 2.1s

%The +2 comes from the numbering on the OST ahead of these commands
for i = 1:n
    if i <= audStimP.rampStps
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(audStimP.tStep) ' NaN {} #Shift ' num2str(i) ' of ' num2str(audStimP.rampStps)];
    elseif i == audStimP.rampStps + 1 
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(audStimP.lenPerVallT) ' NaN {} #Hold for the pitch-shift hold period'];
    elseif i <= 2*audStimP.rampStps + 1
        OST_tline{i+p} = [num2str(i+2) ' ELAPSED_TIME ' num2str(audStimP.tStep) ' NaN {} #Shift ' num2str(i) ' of ' num2str(audStimP.rampStps)];    
    elseif i == 2*audStimP.rampStps + 2
        OST_tline{i+p} = [num2str(i+2) ' OST_END NaN NaN {} #End the dang thing'];
    elseif 1 == 2*audStimP.rampStps + 3
        OST_tline{i+p} = ' ';
    elseif i == 2*audStimP.rampStps + 4
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

%The number of changes to f0 + the hold + the last ONE clean-up line
n = 2*audStimP.rampStps + 1 + 1;

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
    if i <= audStimP.rampStps
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(audStimP.route(i)/100) ', 0.0, 0, 0'];
    elseif i == audStimP.rampStps + 1 
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(audStimP.route(end)/100) ', 0.0, 0, 0'];
    elseif i <= 2*audStimP.rampStps + 1
        PCF_tline{i+p} = [num2str(i+2) ', ' num2str(audStimP.route((2*audStimP.rampStps + 1) - (i-1))/100) ', 0.0, 0, 0'];
    elseif i == 2*audStimP.rampStps + 2
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
axis([0 4 -101 1]); box off;

set(gca, 'FontSize', 16,...
         'FontWeight', 'bold');

end