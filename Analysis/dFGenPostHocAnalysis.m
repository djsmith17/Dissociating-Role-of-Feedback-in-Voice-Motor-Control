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
addressQuest4(dirs, respPer_SomVF, respPer_SomMN, respPer_AudLs)

% Question 5 %%%
addressQuest5(dirs, StatTableAud, StatTableSomMN, StatTableJND)

% Question 6 %%%
addressQuest6(dirs, StatTableJND, StatTableSomVF, StatTableSomMN)

% Question 7 %%%
addressQuest7(dirs, StatTableJND, StatTableAud)

% Question E1 %%%
addressQuestE1(dirs, StatTableJND, StatTableSomMN)

% Question E2 %%%
addressQuestE2(dirs, StatTableJND, StatTableSomMN)
end

function addressQuest4(dirs, respPer_SomVF, respPer_SomMN, respPer_Aud)
% q4: Do participants show similar compensatory respones when only Auditory
% feedback is perturbed? 
plotpos = [10 100];
plotdim = [1200 500];
q2Fig = figure('Color', [1 1 1]);
set(q2Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

fontN = 'Arial';
axisLSize = 12;

color3 = [44 162 95]/255;

plot(respPer_SomVF, 'bo-', 'MarkerFaceColor', 'b'); hold on
plot(respPer_SomMN, 'ro-', 'MarkerFaceColor', 'r'); 
plot(respPer_Aud, 'o-', 'Color', color3, 'MarkerFaceColor', color3);
xlabel('Participant')
ylabel('RespPer (%)')
title('Comparison of Response Percentage between Experimental Conditions')
box off

set(gca,'FontName', fontN,...
    'FontSize', axisLSize,...
    'FontWeight','bold')

legend({'Somato Feedback Pert (No Masking Noise)', 'Somato Feedback Pert (Masking Noise)', 'Aud Feedback Pert'},...
        'Box', 'off',...
        'Edgecolor', [1 1 1],...
        'FontSize', 12,...
        'FontWeight', 'bold')

dirs.quest2FigFile = fullfile(dirs.SavResultsDir, 'Question4.jpg');
export_fig(dirs.quest2FigFile)
end

function addressQuest5(dirs, StatTableAud, StatTableSomMN, StatTableJND)
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
f0            = cell2mat(StatTableJNDLs.f0);

% Perform the correlation
q5AllResponse = [RespPer_Aud RespPer_SomMN];
[corrR, corrP] = corrcoef(q5AllResponse);
q5Sentence = sprintf('Weak positive correlation between Som RespPer and Aud RespPer (n = %d)', length(RespPer_Aud));

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
scatStr.xLabel = 'Aud RespPer (%)';
scatStr.yLabel = 'SomMN RespPer (%)';
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

function addressQuest6(dirs, StatTableJND, StatTableSomVF, StatTableSomMN)
% q6: Is there a relationship between auditory f0 acuity and contribution
% of auditory feedback in the laryngeal perturbation task?

% Currently Expecting Fewer 'Observations' from SomVF and SomMN
I = ismember(StatTableJND.SubjID, StatTableSomVF.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

% I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
% StatTableSomMN(I18,:) = [];
% StatTableSomVF(I18,:) = [];
% StatTableJNDLs(I18,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_Diff  = respPer_SomMN - respPer_SomVF;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = cell2mat(StatTableJNDLs.f0);

q6AllResponse  = [JNDScore respPer_Diff];
[corrR, corrP] = corrcoef(q6AllResponse);
q6Sentence = sprintf('Weak positive correlation between RespPer Diff and JND, (n = %d)', length(respPer_Diff));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(q6AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_Diff;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'RespPer Diff(%) (M-NM)';
scatStr.title  = sprintf('Relationship between Effect of Masking on Response to\nLaryngeal Displacement and Auditory Acuity');
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

function addressQuest7(dirs, StatTableJND, StatTableAud)
% q7: Is there a relationship between the auditory f0 acuity and their
% response to auditory feedback perturbation?

% I18 = strcmp(StatTableAud.SubjID, 'DRF18');
% StatTableAud(I18,:) = [];
% StatTableJND(I18,:) = [];

respPer_Aud = StatTableAud.RespPer;
JNDScore    = StatTableJND.JNDScoreMean;
f0          = cell2mat(StatTableJND.f0);

% Perform the correlation
q7AllResponse = [JNDScore respPer_Aud];
[corrR, corrP] = corrcoef(q7AllResponse);
q7Sentence = sprintf('Weak negative correlation between RespPer and JND (n = %d)', length(JNDScore));

% Perform the correltion controlling for f0
[corrR2, corrP2] = partialcorr(q7AllResponse, f0);

columnNames = {'rho'; 'p_Value'};
rowNames = {'Corr'; 'Partial Corr (f0)'};
allRhos = [round(corrR(1,2), 3); round(corrR2(1,2), 3)];
allPs   = [round(corrP(1,2), 4); round(corrP2(1,2), 4)];
T = table(allRhos, allPs, 'RowNames', rowNames, 'VariableNames', columnNames);

% Organize the output
scatStr.x      = JNDScore;
scatStr.y      = respPer_Aud;
scatStr.xLabel = 'JND Score (cents)';
scatStr.yLabel = 'RespPer (%)';
scatStr.title  = 'Relationship between Response to Auditory Pitch-Shift and and Auditory Acuity';
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

function addressQuestE1(dirs, StatTableJND, StatTableSomMN)
% qE1: Relationship between Response to Laryngeal Displacement and Auditory
% Acuity

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
qE1AllResponse = [JNDScore respPer_SomMN];
[corrR, corrP] = corrcoef(qE1AllResponse);
qE1Sentence = sprintf('Weak positive correlation between RespPer and JND (n = %d)', length(JNDScore));

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
scatStr.yLabel = 'RespPer (%)';
scatStr.title  = 'Relationship between Response to Laryngeal Displacement and Auditory Acuity';
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

function addressQuestE2(dirs, StatTableJND, StatTableSomMN)

% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableJND.SubjID, StatTableSomMN.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

% I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
% StatTableSomMN(I18,:) = [];
% StatTableJNDLs(I18,:) = [];

stimMag_SomMN = StatTableSomMN.StimMag;
JNDScore      = StatTableJNDLs.JNDScoreMean;
f0            = cell2mat(StatTableJNDLs.f0);

% Perform the correlation
qE2AllResponse = [JNDScore stimMag_SomMN];
[corrR, corrP] = corrcoef(qE2AllResponse);
qE2Sentence = sprintf('Weak positive correlation between StimMag and JND (n = %d)', length(JNDScore));

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
scatStr.yLabel = 'StimMag (cents)';
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