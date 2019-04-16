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

dirs.SomResultsFile = fullfile(dirs.SomResultsDir, 'DRF_SomResultsDRF.mat');
dirs.AudResultsFile = fullfile(dirs.AudResultsDir, 'DRF_AudResultsDRF.mat');
dirs.JNDResultsFile = fullfile(dirs.JNDResultsDir, 'DRF_JNDResultsDRF.mat');

load(dirs.SomResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableSom = allSubjRes.statTable;
load(dirs.AudResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableAud = allSubjRes.statTable;
load(dirs.JNDResultsFile) % Returns 'allSubjRes' 'pooledRunStr'
StatTableJND = allSubjRes.statTable;

somVF = strcmp(StatTableSom.AudFB, 'Voice Feedback');
somMN = strcmp(StatTableSom.AudFB, 'Masking Noise');

StatTableSomVF = StatTableSom(somVF, :);
StatTableSomMN = StatTableSom(somMN, :);

% Question 2 %%%
% Currently Expecting Fewer 'Observations' from SomVF and SomMN
I = ismember(StatTableAud.SubjID, StatTableSomVF.SubjID) == 0;
StatTableAudLs = StatTableAud;
StatTableAudLs(I,:) = [];
respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_AudLs = StatTableAudLs.RespPer;

% Ascending Order Compared against SomVF
[~, I] = sort(respPer_SomVF);
respPer_SomVF = respPer_SomVF(I);
respPer_SomMN = respPer_SomMN(I);
respPer_AudLs = respPer_AudLs(I);

% Draw the progression
addressQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_AudLs)

% Question 3 %%%
addressQuest3(dirs, StatTableJND, StatTableSomMN)

% Question 4 %%%
addressQuest4(dirs, StatTableJND, StatTableAud)

% Question 5 %%%
addressQuest5(dirs, StatTableJND, StatTableSomVF, StatTableSomMN)
end

function addressQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_Aud)

plotpos = [10 100];
plotdim = [1200 500];
q2Fig = figure('Color', [1 1 1]);
set(q2Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(respPer_SomVF, 'bo-'); hold on
plot(respPer_SomMN, 'ro-'); 
plot(respPer_Aud, 'go-');
xlabel('Participant')
ylabel('RespPer (%)')
title('Comparison of Response Percentage between Experimental Conditions')
box off

legend({'Somato Feedback Pert (No Masking Noise)', 'Somato Feedback Pert (Masking Noise)', 'Aud Feedback Pert'},...
        'Box', 'off',...
        'Edgecolor', [1 1 1],...
        'FontSize', 12,...
        'FontWeight', 'bold')

dirs.quest2FigFile = fullfile(dirs.SavResultsDir, 'Question2.jpg');
export_fig(dirs.quest2FigFile)
end

function addressQuest3(dirs, StatTableJND, StatTableSomMN)

% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableJND.SubjID, StatTableSomMN.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

% I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
% StatTableSomMN(I18,:) = [];
% StatTableJNDLs(I18,:) = [];

respPer_SomMN = StatTableSomMN.RespPer;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = cell2mat(StatTableJNDLs.f0);

% Perform the correlation
q3AllResponse = [JNDScore respPer_SomMN];
[corrR, corrP] = corrcoef(q3AllResponse);
q3Sentence = sprintf('Weak positive correlation between RespPer and JND Score, R = %.2f, n = %d, P = %.2f', corrR(1,2), length(JNDScore), corrP(1,2));

% Perform the correltion controlling for f0
[corrR, corrP] = partialcorr(q3AllResponse, f0);
q3Sentence2 = sprintf('Weak positive partial correlation between RespPer and JND Score \nwhile controlling for f0, R = %.2f, n = %d, P = %.2f', corrR(1,2), length(JNDScore), corrP(1,2));

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_SomMN;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'RespPer (%)';
scatStr.title  = 'Relationship between Response to Laryngeal Displacement and and Auditory Acuity';
scatStr.sent   = q3Sentence;
scatStr.sent2  = q3Sentence2;
scatStr.qNum   = 3;
scatStr.minY   = min(scatStr.y) - 20;
scatStr.maxY   = max(scatStr.y) + 20;

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuest4(dirs, StatTableJND, StatTableAud)

respPer_Aud = StatTableAud.RespPer;
JNDScore    = StatTableJND.JNDScoreMean;
f0          = cell2mat(StatTableJND.f0);

% Perform the correlation
q4AllResponse = [JNDScore respPer_Aud];
[corrR, corrP] = corrcoef(q4AllResponse);
q4Sentence = sprintf('Weak negative correlation between RespPer and JND Score, R = %.2f, n = %d, P = %.2f', corrR(1,2), length(JNDScore), corrP(1,2));

% Perform the correltion controlling for f0
[corrR, corrP] = partialcorr(q4AllResponse, f0);
q4Sentence2 = sprintf('Weak negative partial correlation between RespPer and JND Score \nwhile controlling for f0, R = %.2f, n = %d, P = %.2f', corrR(1,2), length(JNDScore), corrP(1,2));

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_Aud;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'RespPer (%)';
scatStr.title  = 'Relationship between Response to Auditory Pitch-Shift and and Auditory Acuity';
scatStr.sent   = q4Sentence;
scatStr.sent2  = q4Sentence2;
scatStr.qNum   = 4;
scatStr.minY   = min(scatStr.y) - 20;
scatStr.maxY   = max(scatStr.y) + 20;

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function addressQuest5(dirs, StatTableJND, StatTableSomVF, StatTableSomMN)
% Currently Expecting Fewer 'Observations' from SomVF and SomMN
I = ismember(StatTableJND.SubjID, StatTableSomVF.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_Diff  = respPer_SomMN - respPer_SomVF;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = cell2mat(StatTableJNDLs.f0);

q11AllResponse = [JNDScore respPer_Diff];
[corrR, corrP] = corrcoef(q11AllResponse);
q5Sentence = sprintf('Weak positive correlation between RespPer Diff and JND, R = %.2f, n = %d, P = %.4f', corrR(1,2), length(respPer_Diff), corrP(1,2));

scatStr.x      = JNDScore;
scatStr.y      = respPer_Diff;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'RespPer Diff(%) (M-nM)';
scatStr.title  = sprintf('Relationship between Effect of Masking on Response to\nLaryngeal Displacement and Auditory Acuity');
scatStr.sent   = q5Sentence;
scatStr.qNum   = 5;
scatStr.minY   = min(scatStr.y) - 20;
scatStr.maxY   = max(scatStr.y) + 20;

% Draw the scatter plot
drawScatterCorr(dirs, scatStr)
end

function drawScatterCorr(dirs, scatStr)
% drawScatterCorr(dirs, respPer_Aud, JNDScore, sentence) draws a
% scatterplot revealing the relationship between two independent variables.
% This also shares the R statistic of the relationship, and the p-value.

plotpos = [10 100];
plotdim = [900 500];
scatFig = figure('Color', [1 1 1]);
set(scatFig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

fontN = 'Arial';
axisLSize = 14;

plot(scatStr.x, scatStr.y, 'ko', 'MarkerSize', 10)
xlabel(scatStr.xLabel)
ylabel(scatStr.yLabel)
title(scatStr.title)
axis([0 70 scatStr.minY scatStr.maxY])
box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

annotation('Textbox', [0.15 0.825 0.1 0.1], ...
           'String', scatStr.sent,...
           'FontName', fontN,...
           'FontSize', 12,...
           'FontWeight','bold',...
           'EdgeColor', 'none')
       
annotation('Textbox', [0.15 0.785 0.1 0.1], ...
           'String', scatStr.sent2,...
           'FontName', fontN,...
           'FontSize', 12,...
           'FontWeight','bold',...
           'EdgeColor', 'none')

dirs.scatFigFile = fullfile(dirs.SavResultsDir, ['Question' num2str(scatStr.qNum) '.jpg']);
export_fig(dirs.scatFigFile)
end