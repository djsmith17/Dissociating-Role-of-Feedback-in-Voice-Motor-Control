clear

resultsFolder = 'C:\Users\djsmith\Documents\MATLAB\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Results\Pilot28\LDDdiag4';

run = 'LDDDiag4';

AudFile = fullfile(resultsFolder, ['Pilot28' run 'ResultsAudDRF.mat']);
load(AudFile)
auRes = res;

PraatFile = fullfile(resultsFolder, ['Pilot28' run 'ResultsPraatDRF.mat']);
load(PraatFile)
prRes = res;

auOn    = auRes.audioMf0MeanPert(:,1);
auOnCI  = auRes.audioMf0MeanPert(:,2);
auOf    = auRes.audioMf0MeanPert(:,3);
auOfCI  = auRes.audioMf0MeanPert(:,4);
auNumTrial = auRes.numPertTrialsPP;
auetMH = auRes.etMH;

prOn    = prRes.audioMf0MeanPert(:,1);
prOnCI  = prRes.audioMf0MeanPert(:,2);
prOf    = prRes.audioMf0MeanPert(:,3);
prOfCI  = prRes.audioMf0MeanPert(:,4);
prNumTrial = prRes.numPertTrialsPP;
pretMH = prRes.etMH;

auLen = length(auOn);
prLen = length(prOn);

if auLen ~= prLen
    auOn    = resample(auOn, prLen, auLen);
    auOnCI  = resample(auOnCI, prLen, auLen);
    auOf    = resample(auOf, prLen, auLen);
    auOfCI  = resample(auOfCI, prLen, auLen);
end

errOn   = prOn - auOn;
errOnCI = prOnCI - auOnCI;
errOf   = prOf - auOf;
errOfCI = prOfCI - auOfCI;

SQEOn   = errOn.^2;
SQEOnCI = errOnCI.^2;
SQEOf   = errOf.^2;
SQEOfCI = errOfCI.^2;

MSEOn    = mean(SQEOn);
MSEOnCI  = mean(SQEOnCI);
MSEOf    = mean(SQEOf);
MSEOfCI  = mean(SQEOfCI);

RMSEOn   = sqrt(MSEOn);
RMSEOnCI = sqrt(MSEOnCI);
RMSEOf   = sqrt(MSEOf);
RMSEOfCI = sqrt(MSEOfCI);

fprintf('%f %f %f %f\n', RMSEOn, RMSEOnCI, RMSEOf, RMSEOfCI)

figure
plot(auOn)
hold on
plot(prOn,'r')
legend(['Audapter: ' num2str(auNumTrial) ' trials, ' num2str(auetMH) ' min'],...
       ['Praat: ' num2str(prNumTrial) ' trials, ' num2str(pretMH) ' min'])