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

respPer_SomVF = StatTableSom.RespPer(somVF);
respPer_SomMN = StatTableSom.RespPer(somMN);
respPer_Aud   = StatTableAud.RespPer;
JNDScore      = StatTableJND.JNDScoreMean;

% Ascending Order Compared against SomVF
[~, I] = sort(respPer_SomVF);
respPer_SomVF = respPer_SomVF(I);
respPer_SomMN = respPer_SomMN(I);
respPer_Aud   = respPer_Aud(I);
JNDScore      = JNDScore(I);

allRespPer = [respPer_SomMN respPer_Aud JNDScore];
[corrR, corrP] = corrcoef(allRespPer);

q3Sentence = sprintf('Weak positive correlation between RespPer and JND Score, R = %.2f, P = %.2f', corrR(1,3), corrP(1,3));
q4Sentence = sprintf('Weak negative correlation between RespPer and JND Score, R = %.2f, P = %.2f', corrR(2,3), corrP(2,3));

drawQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_Aud)
drawQuest3(dirs, respPer_SomMN, JNDScore, q3Sentence)
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