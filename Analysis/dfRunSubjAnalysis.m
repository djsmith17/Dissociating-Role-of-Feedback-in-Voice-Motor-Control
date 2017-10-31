function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
%Plot Toggles. This could eventually become an input variable
PltTgl.ForceSensor     = 0; %Voltage trace of force sensor signal
PltTgl.IntraTrial_T    = 0; %SPL trace of individual trial
PltTgl.IntraTrial_f0   = 0; %f0 trace of individual trial
PltTgl.InterTrial_f0   = 1; %Average f0 trace over all trials of a run
PltTgl.InterRun_f0       = 1; %Average f0 trace over all runs analyzed
PltTgl.InterTrial_AudRes = 1; %Average f0 response trace to auditory pert trials of a run
PltTgl.InterRun_AudRes   = 1; %Average f0 response trace to auditory pert over all runs analyzed
PltTgl.InterTrial_Force  = 0;
PltTgl.InterRun_Force    = 0;
PltTgl.svInflaRespRoute  = 0;

AVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participant  = 'Pilot0'; %List of multiple participants.
AVar.run          = 'SF3';

dirs = dfDirs(AVar.project);

dirs.SavFileDir    = fullfile(dirs.SavData, AVar.participant, AVar.run, [AVar.participant AVar.run 'DRF.mat']); %Where to find data
dirs.SavResultsDir = fullfile(dirs.Results, AVar.participant, AVar.run); %Where to save results
dirs.InflaRespFile = fullfile(dirs.SavData, AVar.participant, [AVar.participant '_AveInflaResp.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

fprintf('Loading Files for %s %s\n', AVar.participant, AVar.run)
load(dirs.SavFileDir)

[auAn, auRes] = dfAnalysisAudapter(DRF.expParam, DRF.rawData, DRF.DAQin);
[niAn, niRes] = dfAnalysisNIDAQ(DRF.expParam, DRF.DAQin);

fprintf('Saving Results for %s %s\n', AVar.participant, AVar.run)
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [AVar.participant AVar.run 'ResultsDRF.mat']);
save(dirs.SavResultsFile, 'auAn', 'auRes', 'niAn', 'niRes')
% 
% drawDAQAll(niAn, 2, dirs.SavResultsDir, 1)
% drawInterTrialf0(auRes.timeSec, auRes.meanTrialf0_St, auRes.meanTrialf0_Sp, auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
% drawAllTrialf0(auRes.time, auRes.allTrialf0, auRes.runTrialOrder, auAn.trigsT, auRes.f0Limits, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
% drawAudResp_AllTrial(auRes, auAn.curSess, DRF.expParam.curSess, dirs.SavResultsDir)
% 
% drawAudResp_InterTrial(auRes.timeSec, auRes.meanTrialf0_St(:,:,2), auRes.meanTrialf0_Sp(:,:,2), auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
%          
% if PltTgl.svInflaRespRoute == 1
%     InflaRespRoute = CalcInflationResponse(auAn, auRes.meanTrialf0b, auRes.meanTrialf0_St, auRes.InflaRespLimits, dirs.SavResultsDir);
%     tStep = auAn.tStep;
%     save(dirs.InflaRespFile, 'InflaRespRoute', 'tStep')
% end
end

function InflaRespRoute = CalcInflationResponse(auAn, meanTrialf0b, meanRunsf0, limits, plotFolder)
%This calculates the shape of the change in f0 in response to the
%perturbation onset

time = auAn.time;
curExp = auAn.curExp;
tStep  = auAn.tStep;

disp('Calculating Inflation Response Route')

meanPertMicf0 = meanRunsf0(:,1,2);           %Grabbing the mean mic f0 for all perturbed trials
postOnset     = find(time > 0.5); %Trials are centered at 0.5s before inflation. 
[~, ind]      = min(meanPertMicf0);          %Find the lowest point in whole mean signal

timeFram = postOnset(1):ind;
InflaRespRoute = meanPertMicf0(timeFram);

plotpos = [200 400];
plotdim = [600 600];
InflaRespFig = figure('Color',[1 1 1]);
set(InflaRespFig, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

t = 0:tStep:(tStep*(length(timeFram)-1)); %meanPertMicf0(timeFram,1) - meanPertMicf0(postOnset(1));
plot(t, InflaRespRoute, 'blue', 'LineWidth',3)

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({'Average Acoustic Response to Inflation'; [curExp '   f0: ' num2str(meanTrialf0b) 'Hz']},...
                     'FontSize', 18,...
                     'FontWeight', 'bold')
axis(limits); box off

set(gca, 'FontSize', 16,...
         'FontWeight','bold');

plTitle = [curExp '_Inflation Response Route.jpg'];
saveFileName = fullfile(plotFolder, plTitle);
export_fig(saveFileName)
end