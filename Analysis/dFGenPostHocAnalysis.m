function dfGenPostHocAnalysis()

close all

pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.phExp         = 'DRF_PostHoc';

pA.SomExp     = 'DRF_Som';
pA.AudExp     = 'DRF_Aud';
pA.JNDExp     = 'DRF_JND';

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
drawQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_AudLs)

% Question 2.1
respPerDiff = respPer_SomVF - respPer_SomMN;
q21AllResponse = [respPerDiff respPer_SomMN];
[corrR, corrP] = corrcoef(q21AllResponse);
q21Sentence = sprintf('Strong negative correlation between RespPer and JND Score, R = %.2f, P = %.4f', corrR(1,2), corrP(1,2));

% Draw the scatter plot
drawQuest21(dirs, respPer_SomMN, respPerDiff, q21Sentence)

% Question 3 %%%
% Currently Expecting Fewer 'Observations' from SomMN
I = ismember(StatTableJND.SubjID, StatTableSomMN.SubjID) == 0;
StatTableJNDLs = StatTableJND;
StatTableJNDLs(I,:) = [];

% I18 = strcmp(StatTableSomMN.SubjID, 'DRF18');
% StatTableSomMN(I18,:) = [];
% StatTableJNDLs(I18,:) = [];

respPer_SomMN = StatTableSomMN.RespPer;
JNDScoreLs    = StatTableJNDLs.JNDScoreMean;

% Perform the correlation
q3AllResponse = [JNDScoreLs respPer_SomMN];
[corrR, corrP] = corrcoef(q3AllResponse);
q3Sentence = sprintf('Weak positive correlation between RespPer and JND Score, R = %.2f, P = %.2f', corrR(1,2), corrP(1,2));

% Draw the scatter plot
drawQuest3(dirs, respPer_SomMN, JNDScoreLs, q3Sentence)

% Question 4 %%%

% I18 = strcmp(StatTableAud.SubjID, 'DRF18');
% StatTableAud(I18,:) = [];
% StatTableJND(I18,:) = [];

respPer_Aud = StatTableAud.RespPer;
JNDScore    = StatTableJND.JNDScoreMean;

% Perform the correlation
q4AllResponse = [JNDScore respPer_Aud];
[corrR, corrP] = corrcoef(q4AllResponse);
q4Sentence = sprintf('Weak negative correlation between RespPer and JND Score, R = %.2f, P = %.2f', corrR(1,2), corrP(1,2));

% Draw the scatter plot
drawQuest4(dirs, respPer_Aud, JNDScore, q4Sentence)
end

function drawQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_Aud)

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

legend({'Somato Feedback Pert (Voice Feedback)', 'Somato Feedback Pert (Masking Noise)', 'Aud Feedback Pert'},...
        'Box', 'off',...
        'Edgecolor', [1 1 1],...
        'FontSize', 12,...
        'FontWeight', 'bold')

dirs.quest2FigFile = fullfile(dirs.SavResultsDir, 'Question2.jpg');
export_fig(dirs.quest2FigFile)
end

function drawQuest21(dirs, respPer_Som, JNDScore, sentence)

plotpos = [10 100];
plotdim = [900 500];
q21Fig = figure('Color', [1 1 1]);
set(q21Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(JNDScore, respPer_Som, 'ko')
xlabel('Diff RespPer (nMN-MN)')
ylabel('RespPer (%)')
title('Relationship between Somato Pert RespPer (Masking Noise) and Effect of Masking')
box off

annotation('Textbox', [0.15 0.84 0.1 0.1], ...
           'String', sentence,...
           'EdgeColor', 'none')

dirs.quest21FigFile = fullfile(dirs.SavResultsDir, 'Question21.jpg');
export_fig(dirs.quest21FigFile)
end

function drawQuest3(dirs, respPer_Som, JNDScore, sentence)

plotpos = [10 100];
plotdim = [900 500];
q3Fig = figure('Color', [1 1 1]);
set(q3Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(JNDScore, respPer_Som, 'ko')
xlabel('JND Score')
ylabel('RespPer (%)')
title('Relationship between Somato Pert RespPer and JND Score')
box off

annotation('Textbox', [0.15 0.84 0.1 0.1], ...
           'String', sentence,...
           'EdgeColor', 'none')

dirs.quest3FigFile = fullfile(dirs.SavResultsDir, 'Question3.jpg');
export_fig(dirs.quest3FigFile)
end

function drawQuest4(dirs, respPer_Aud, JNDScore, sentence)

plotpos = [10 100];
plotdim = [900 500];
q4Fig = figure('Color', [1 1 1]);
set(q4Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(JNDScore, respPer_Aud, 'ko')
xlabel('JND Score')
ylabel('RespPer (%)')
title('Relationship between Audio Pert RespPer and JND Score')
box off

annotation('Textbox', [0.15 0.84 0.1 0.1], ...
           'String', sentence,...
           'EdgeColor', 'none')

dirs.quest4FigFile = fullfile(dirs.SavResultsDir, 'Question4.jpg');
export_fig(dirs.quest4FigFile)
end