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

load(dirs.SomResultsFile) %Returns allSubjRes pooledRunStr
SomStatTable = allSubjRes.statTable;
load(dirs.AudResultsFile) %Returns allSubjRes pooledRunStr
AudStatTable = allSubjRes.statTable;
load(dirs.JNDResultsFile) %Returns allSubjRes pooledRunStr
JNDStatTable = allSubjRes.statTable;

somVF = strcmp(SomStatTable.AudFB, 'Voice Feedback');
somMN = strcmp(SomStatTable.AudFB, 'Masking Noise');

respPer_SomVF = SomStatTable.RespPer(somVF);
respPer_SomMN = SomStatTable.RespPer(somMN);
respPer_Aud   = AudStatTable.RespPer;
JNDScore      = JNDStatTable.JNDScoreMean;

[~, I] = sort(respPer_SomVF);
respPer_SomVF = respPer_SomVF(I);
respPer_SomMN = respPer_SomMN(I);
respPer_Aud   = respPer_Aud(I);
JNDScore      = JNDScore(I);

drawQuest2(dirs, respPer_SomVF, respPer_SomMN, respPer_Aud)
drawQuest3(dirs, respPer_SomVF, JNDScore)
drawQuest4(dirs, respPer_Aud, JNDScore)
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

function drawQuest3(dirs, respPer_SomVF, JNDScore)

plotpos = [10 100];
plotdim = [900 500];
q3Fig = figure('Color', [1 1 1]);
set(q3Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(JNDScore, respPer_SomVF, 'ko')
xlabel('JND Score')
ylabel('RespPer (%)')
title('Relationship between Somato Pert RespPer and JND Score')

dirs.quest3FigFile = fullfile(dirs.SavResultsDir, 'Question3.jpg');
export_fig(dirs.quest3FigFile)
end

function drawQuest4(dirs, respPer_Aud, JNDScore)

plotpos = [10 100];
plotdim = [900 500];
q4Fig = figure('Color', [1 1 1]);
set(q4Fig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(JNDScore, respPer_Aud, 'ko')
xlabel('JND Score')
ylabel('RespPer (%)')
title('Relationship between Audio Pert RespPer and JND Score')

dirs.quest4FigFile = fullfile(dirs.SavResultsDir, 'Question4.jpg');
export_fig(dirs.quest4FigFile)
end