
close all; clear all;

%Plotting Variables
plotpos    = [100 100];
plotdim    = [600 800];
DIVAResp   = figure('Color', [1 1 1]);
set(DIVAResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

fontN        = 'Arial';
legAnnoFSize = 25;
titleFSize   = 20;
axisLSize    = 30;
lineThick    = 4;

fs = 10000; %samp/s
t = 0:1/fs:4-1/fs;

tOnset  = 1*fs;
tOffset = 3*fs;

tPitchRespOnset = (0:.50*fs) + tOnset;
lenPitchResp = length(tPitchRespOnset);

fpitch = 1/(2*lenPitchResp);

pertSig = zeros(size(t));
pertSig(tOnset:tOffset) = 1;
pertSig = pertSig + 3;

pitchSig = zeros(size(t));
pitchSig(tPitchRespOnset) = -sin(2*pi*fpitch*(1:lenPitchResp));
pitchSig = pitchSig - 2;

plot(t, pertSig, 'k', 'LineWidth', lineThick)
hold on
plot(t, pitchSig, 'r', 'LineWidth', lineThick)

title('Laryngeal Displacement Dynamics', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis([0 4 -5 5]); box off;
% axis off