function dfGenPostHocAnalysis()
% This function sets up the relationships between participant responses in
% my dissertation project. Specifically, this looks at the responses from
% DRF_Som, DRF_Aud, & DRF_JND
close all

pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.phExp         = 'DRF_PostHoc';

pA.SomExp     = 'DRF_Som'; % Somatosensory Feedback Pert Pooled Analysis
pA.AudExp     = 'DRF_Aud'; % Auditory Feedback Pert Pooled Analysis
pA.JNDExp     = 'DRF_JND'; % f0 Acuity Pooled Analysis

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.phExp);
dirs.SomResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.SomExp);
dirs.AudResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.AudExp);
dirs.JNDResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.JNDExp);

% Write Results to Text File
textFileName = fullfile(dirs.SavResultsDir, 'PostHocStats.txt');
fid = fopen(textFileName, 'w');

dirs.SomResultsFile = fullfile(dirs.SomResultsDir, 'DRF_SomResultsDRF.mat');
dirs.AudResultsFile = fullfile(dirs.AudResultsDir, 'DRF_AudResultsDRF.mat');
dirs.JNDResultsFile = fullfile(dirs.JNDResultsDir, 'DRF_JNDResultsDRF.mat');

load(dirs.SomResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableSom       = allSubjRes.statTable;
StatTableSomSingle = allSubjRes.statTableSingle;
pooledRunStrSom    = pooledRunStr;
load(dirs.AudResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableAud       = allSubjRes.statTable;
StatTableAudSingle = allSubjRes.statTableSingle;
load(dirs.JNDResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableJND = allSubjRes.statTable;

somVF = strcmp(StatTableSom.AudFB, 'Voice Feedback');
somMN = strcmp(StatTableSom.AudFB, 'Masking Noise');

StatTableSomVF = StatTableSom(somVF, :);
StatTableSomMN = StatTableSom(somMN, :);

StatTableSinSomVF = StatTableSomSingle(somVF, :);
StatTableSinSomMN = StatTableSomSingle(somMN, :);

% Question 4 %%%
StatsOrg_DRF_Som_Aud_RMANOVA_covar(dirs, StatTableSomVF, StatTableSomMN, StatTableAud, StatTableJND)

% Question 5 %%%
addressQuest5(dirs, fid, StatTableAud, StatTableSomMN, StatTableJND)

% Question 6 %%%
addressQuest6(dirs, fid, StatTableJND, StatTableSomVF, StatTableSomMN)

% Question 7 %%%
addressQuest7(dirs, fid, StatTableJND, StatTableAud)

% Question E1 %%%
addressQuestE1(dirs, fid, StatTableJND, StatTableSomMN)

% Question E2 %%%
addressQuestE2(dirs, fid, StatTableJND, StatTableSomMN)

% Question E3 %%%
addressQuestE3(dirs, fid, StatTableJND)

% Question E4 and E5 %%%
addressQuestE4(dirs, fid, pooledRunStrSom)
addressQuestE5(dirs, fid, pooledRunStrSom)

fclose(fid);
end

function addressQuest5(dirs, fid, StatTableAud, StatTableSomMN, StatTableJND)
% q5: Is there a relationship between individual response magnitudes to
% somatosensory and auditory perturbations?

% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableAud.SubjID, StatTableSomMN.SubjID) == 0;
StatTableAudLs = StatTableAud;
StatTableAudLs(I,:) = [];
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

% I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
% StatTableSomMN(I18,:) = [];
% StatTableJNDLs(I18,:) = [];

RespPer_SomMN = StatTableSomMN.RespPer;
RespPer_Aud   = StatTableAudLs.RespPer;
f0            = StatTableJNDLs.f0;

% Perform the correlation
q5AllResponse = [RespPer_Aud RespPer_SomMN];
[corrR, corrP] = corrcoef(q5AllResponse);
q5Sentence = sprintf('Weak positive correlation between Som RespPer and Aud RespPer (n = %d)', length(RespPer_Aud));

fprintf(fid, 'Relationship between responses to Somatosensory and Auditory Perturbations\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n\n', corrR(1,2), corrP(1,2), length(RespPer_Aud));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(q5AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = RespPer_Aud;
scatStr.y      = RespPer_SomMN;
scatStr.xLabel = 'Audio. Response Percentage (%)';
scatStr.yLabel = 'Somato. Response Percentage (%)';
scatStr.title  = 'Relationship between Responses to Somatosensory and Auditory Perts';
scatStr.sent   = q5Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 5;
scatStr.minX   = min(scatStr.x) - 10;
scatStr.maxX   = max(scatStr.x) + 10;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.38 0.23];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuest6(dirs, fid, StatTableJND, StatTableSomVF, StatTableSomMN)
% q6: Is there a relationship between auditory f0 acuity and contribution
% of auditory feedback in the laryngeal perturbation task?

% Currently Expecting Fewer 'Observations' from SomVF and SomMN
I = ismember(StatTableJND.SubjID, StatTableSomVF.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_Diff  = respPer_SomVF - respPer_SomMN; % Not Masked - Masked
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = StatTableJNDLs.f0;
RMS           = StatTableJNDLs.audioRMS;

q6AllResponse  = [JNDScore respPer_Diff];
[corrR, corrP] = corrcoef(q6AllResponse);
q6Sentence = sprintf('Weak positive correlation between RespPer Diff and JND, (n = %d)', length(respPer_Diff));

fprintf(fid, 'Relationship between Effect of Masking on Response to Laryngeal Perturbations and Auditory Acuity\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n', corrR(1,2), corrP(1,2), length(JNDScore));

%Outlier
I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
StatTableSomMN(I18,:) = [];
StatTableSomVF(I18,:) = [];
StatTableJNDLs(I18,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_Diff_Out  = respPer_SomVF - respPer_SomMN; % Not Masked - Masked
RMS_Out           = StatTableJNDLs.audioRMS;
JNDScore_Out      = StatTableJNDLs.JNDScoreMean;
q6AllResponse_Out = [JNDScore_Out respPer_Diff_Out];
[corrR, corrP] = corrcoef(q6AllResponse_Out);

fprintf(fid, 'With the outlier removed\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n\n', corrR(1,2), corrP(1,2), length(JNDScore_Out));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(q6AllResponse, RMS);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_Diff;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'Response Percentage Diff(%) (Wo - W)';
scatStr.title  = sprintf('Relationship between Effect of Masking on Response to\nLaryngeal Perturbation and Auditory Acuity');
scatStr.sent   = q6Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 6;
scatStr.minX   = 0;
scatStr.maxX   = 70;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.42 0.26];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuest7(dirs, fid, StatTableJND, StatTableAud)
% q7: Is there a relationship between the auditory f0 acuity and their
% response to auditory feedback perturbation?

respPer_Aud = StatTableAud.RespPer;
JNDScore    = StatTableJND.JNDScoreMean;
f0          = StatTableJND.f0;
RMS         = StatTableJND.audioRMS;

% Perform the correlation
q7AllResponse = [JNDScore respPer_Aud];
[corrR, corrP] = corrcoef(q7AllResponse);
q7Sentence = sprintf('Weak negative correlation between RespPer and JND (n = %d)', length(JNDScore));

fprintf(fid, 'Relationship between Response to Auditory Perturbations and Auditory Acuity\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n', corrR(1,2), corrP(1,2), length(JNDScore));

%Outlier
I18 = strcmp(StatTableAud.SubjID, 'DRF18');
StatTableAud(I18,:) = [];
StatTableJND(I18,:) = [];

respPer_Aud_Out = StatTableAud.RespPer;
JNDScore_Out    = StatTableJND.JNDScoreMean;

q7AllResponse_Out = [JNDScore_Out respPer_Aud_Out];
[corrR, corrP] = corrcoef(q7AllResponse_Out);
fprintf(fid, 'With the outlier removed\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n\n', corrR(1,2), corrP(1,2), length(JNDScore_Out));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(q7AllResponse, RMS);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_Aud;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'Response Percentage (%)';
scatStr.title  = 'Relationship between Response to Auditory Feedback Perturbation and Auditory Acuity';
scatStr.sent   = q7Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 7;
scatStr.minX   = 0;
scatStr.maxX   = 70;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.50 0.69];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuestE1(dirs, fid, StatTableJND, StatTableSomMN)
% qE1: Relationship between Response to Laryngeal Displacement and Auditory
% Acuity

% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableJND.SubjID, StatTableSomMN.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

respPer_SomMN = StatTableSomMN.RespPer;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = StatTableJNDLs.f0;
RMS           = StatTableJNDLs.audioRMS;

% Perform the correlation
qE1AllResponse = [JNDScore respPer_SomMN];
[corrR, corrP] = corrcoef(qE1AllResponse);
qE1Sentence = sprintf('Weak positive correlation between RespPer and JND (n = %d)', length(JNDScore));

fprintf(fid, 'Relationship between Response to Laryngeal Perturbations and Auditory Acuity\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n', corrR(1,2), corrP(1,2), length(JNDScore));

%Outlier
I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
StatTableSomMN(I18,:) = [];
StatTableJNDLs(I18,:) = [];

respPer_SomMN_Out = StatTableSomMN.RespPer;
JNDScore_Out      = StatTableJNDLs.JNDScoreMean;

qE1AllResponse_Out = [JNDScore_Out respPer_SomMN_Out];
[corrR, corrP] = corrcoef(qE1AllResponse_Out);
fprintf(fid, 'With the outlier removed\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n\n', corrR(1,2), corrP(1,2), length(JNDScore));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(qE1AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_SomMN;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'Response Percentage (%)';
scatStr.title  = 'Relationship between Response to Laryngeal Perturbations and Auditory Acuity';
scatStr.sent   = qE1Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 11;
scatStr.minX   = 0;
scatStr.maxX   = 70;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.52 0.23];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuestE2(dirs, fid, StatTableJND, StatTableSomMN)

% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableJND.SubjID, StatTableSomMN.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

stimMag_SomMN = StatTableSomMN.StimMag;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = StatTableJNDLs.f0;

% Perform the correlation
qE2AllResponse = [JNDScore stimMag_SomMN];
[corrR, corrP] = corrcoef(qE2AllResponse);
qE2Sentence = sprintf('Weak positive correlation between StimMag and JND (n = %d)', length(JNDScore));

fprintf(fid, 'Relationship between Stimulus Magnitude by Laryngeal Perturbations and Auditory Acuity\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n', corrR(1,2), corrP(1,2), length(JNDScore));

I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
StatTableSomMN(I18,:) = [];
StatTableJNDLs(I18,:) = [];

stimMag_SomMN_Out = StatTableSomMN.StimMag;
JNDScore_Out      = StatTableJNDLs.JNDScoreMean;

qE2AllResponse_Out = [JNDScore_Out stimMag_SomMN_Out];
[corrR, corrP] = corrcoef(qE2AllResponse_Out);

fprintf(fid, 'With the outlier removed\n');
fprintf(fid, 'r = %0.5f, p = %0.5f, n = %d\n\n', corrR(1,2), corrP(1,2), length(JNDScore));


% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(qE2AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = stimMag_SomMN;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'Stimulus Magnitude (cents)';
scatStr.title  = 'Relationship between Magnitude of Laryngeal Perturbation and Auditory Acuity';
scatStr.sent   = qE2Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 22;
scatStr.minX   = 0;
scatStr.maxX   = 70;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.48 0.73];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuestE3(dirs, fid, StatTableJND)

% I18 = strcmp(StatTableJND.SubjID, 'DRF18');
% StatTableJND(I18,:) = [];

JNDScore      = StatTableJND.JNDScoreMean;
f0            = StatTableJND.f0;
RMS           = StatTableJND.audioRMS;

% Perform the correlation
qE3AllResponse = [JNDScore RMS];
[corrR, corrP] = corrcoef(qE3AllResponse);
qE3Sentence = sprintf('Weak positive correlation between StimMag and JND (n = %d)', length(JNDScore));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(qE3AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = RMS;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'Loudness Level (dB)';
scatStr.title  = 'Relationship between JND Score and Token loudness';
scatStr.sent   = qE3Sentence;
scatStr.color  = 'k';
scatStr.Table  = T;
scatStr.qNum   = 33;
scatStr.minX   = 0;
scatStr.maxX   = 70;
scatStr.minY   = min(scatStr.y) - 10;
scatStr.maxY   = max(scatStr.y) + 10;
scatStr.winPos = [0.48 0.73];

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuestE4(dirs, fid, pooledRunStrSom)
% Question E4% Perturbation Lengths

% measures to consider
meas    = {'PertDur'};
numMeas = length(meas);
cond    = {'VF', 'MN'};
numCond = length(cond);

pertDurVF = [];
pertDurMN = [];

for ii = 1:length(pooledRunStrSom)
   pertDurVF = cat(2, pertDurVF, pooledRunStrSom(ii).pertLengths{1});
   pertDurMN = cat(2, pertDurMN, pooledRunStrSom(ii).pertLengths{2}); 
end

pertDurVFMean = mean(pertDurVF);
pertDurMNMean = mean(pertDurMN);

[~, P, ~, STATS] = ttest([pertDurVFMean'-pertDurMNMean']);

pertDurVFMeanMean = mean(pertDurVFMean);
pertDurMNMeanMean = mean(pertDurMNMean);

pertDurVFSTD = std(pertDurVFMean);
pertDurMNSTD = std(pertDurMNMean);

figure
boxplot(pertDurVF)
xlabel('Participant')
ylabel('Perturbation Length (s)')
title({'Without Masking Condition', 'Mean = 1.257s'})

figure
boxplot(pertDurMN)
xlabel('Participant')
ylabel('Perturbation Length (s)')
title({'With Masking Condition', 'Mean = 1.260s'})

fprintf('The mean duration times (With Masking: %.2f +/- %.2f; Without Masking: %.2f +/- %.2f) did not differ significantly by condition (duration: t(%d)= %.2f, p = %.3f).\n', pertDurMNMeanMean,...
                                                                                                                                                                                pertDurMNSTD,...
                                                                                                                                                                                pertDurVFMeanMean,...
                                                                                                                                                                                pertDurVFSTD,...
                                                                                                                                                                                STATS.df,...
                                                                                                                                                                                STATS.tstat,...
                                                                                                                                                                                P)

end

function addressQuestE5(dirs, fid, pooledRunStrSom)
% Question E5% Perturbation Lengths

% measures to consider
meas    = {'PertOnsetTime'};
numMeas = length(meas);
cond    = {'VF', 'MN'};
numCond = length(cond);

pertOnsetTimesVF = [];
pertOnsetTimesMN = [];

for ii = 1:length(pooledRunStrSom)
   pertOnsetTimesVF = cat(2, pertOnsetTimesVF, pooledRunStrSom(ii).pertOnsetTimes{1});
   pertOnsetTimesMN = cat(2, pertOnsetTimesMN, pooledRunStrSom(ii).pertOnsetTimes{2}); 
end

pertOnsetTimesVFMean = mean(pertOnsetTimesVF);
pertOnsetTimesMNMean = mean(pertOnsetTimesMN);

[~, P, ~, STATS] = ttest([pertOnsetTimesVFMean'-pertOnsetTimesMNMean']);

pertOnsetTimesVFMeanMean = mean(pertOnsetTimesVFMean);
pertOnsetTimesMNMeanMean = mean(pertOnsetTimesMNMean);

pertOnsetTimesVFSTD = std(pertOnsetTimesVFMean);
pertOnsetTimesMNSTD = std(pertOnsetTimesMNMean);

figure
boxplot(pertOnsetTimesVF)
xlabel('Participant')
ylabel('Perturbation Length (s)')
title({'Without Masking Condition', 'Mean = 2.10s'})
figure
boxplot(pertOnsetTimesMN)
xlabel('Participant')
ylabel('Perturbation Length (s)')
title({'With Masking Condition', 'Mean = 2.13s'})

fprintf('The mean onset times (With Masking: %.2f +/- %.2f; Without Masking: %.2f +/- %.2f) did not differ significantly by condition (onset: t(%d)= %.2f, p = %.3f).\n', pertOnsetTimesMNMeanMean,...
                                                                                                                                                                          pertOnsetTimesMNSTD,...
                                                                                                                                                                          pertOnsetTimesVFMeanMean,...
                                                                                                                                                                          pertOnsetTimesVFSTD,...
                                                                                                                                                                          STATS.df,...
                                                                                                                                                                          STATS.tstat,...
                                                                                                                                                                          P)
end

function drawScatterCorr(dirs, scatStr)
% drawScatterCorr(dirs, respPer_Aud, JNDScore, sentence) draws a
% scatterplot revealing the relationship between two independent variables.
% This also shares the R statistic of the relationship, and the p-value.

plotpos = [10 100];
plotdim = [900 500];
scatFig = figure('Color', [1 1 1]);
set(scatFig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

fontN = 'Times New Roman';
axisLSize = 14;

plot(scatStr.x, scatStr.y, 'o', 'MarkerSize', 10, 'Color', scatStr.color, 'MarkerFaceColor', scatStr.color)
xlabel(scatStr.xLabel)
ylabel(scatStr.yLabel)
title(scatStr.title)
axis([scatStr.minX scatStr.maxX scatStr.minY scatStr.maxY])
box off

corrSentX  = scatStr.winPos(1);
corrSentY1 = scatStr.winPos(2);
tableH     = 0.127;
tableW     = 0.377;
tableY     = corrSentY1 - tableH + 0.06;
tableX     = corrSentX + 0.04;

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

% annotation('Textbox', [corrSentX corrSentY1 0.1 0.1], ...
%            'String', scatStr.sent,...
%            'FontName', fontN,...
%            'FontSize', 10,...
%            'FontWeight','bold',...
%            'EdgeColor', 'none',...
%            'Color', scatStr.color)
%        
% ScatT = uitable('Data', scatStr.Table{:,:},'ColumnName',scatStr.Table.Properties.VariableNames,...
%         'RowName',scatStr.Table.Properties.RowNames, 'Units', 'Normalized', 'Position',[tableX, tableY, tableW, tableH],...
%         'FontWeight', 'Bold');

dirs.scatFigFile = fullfile(dirs.SavResultsDir, ['Question' num2str(scatStr.qNum) '.jpg']);
export_fig(dirs.scatFigFile, '-nocrop')
end