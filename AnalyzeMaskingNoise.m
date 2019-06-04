clear
close all
[AudMask, fsA] = audioread('C:\GitHub\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Presentation\PrelimFiles\SSN.wav');
[TacMask, fsT] = audioread('C:\GitHub\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Presentation\PrelimFiles\SSN.wav');

N = length(AudMask);
xdft = fft(AudMask);
xdft = xdft(1:N/2+1);
psdxA = (1/(fsA*N)) * abs(xdft).^2;
psdxA(2:end-1) = 2*psdxA(2:end-1);
freqA = 0:fsA/N:fsA/2;

N = length(TacMask);
xdft = fft(TacMask);
xdft = xdft(1:N/2+1);
psdxT = (1/(fsT*N)) * abs(xdft).^2;
psdxT(2:end-1) = 2*psdxT(2:end-1);
freqT = 0:fsT/N:fsT/2;


plotpos = [10 100];
plotdim = [1600 600];

MaskFig = figure('Color', [1 1 1]);
set(MaskFig, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

plot(freqA,10*log10(psdxA), 'r')
hold on
% plot(freqT,10*log10(psdxT), 'k')

grid on
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Comparison of Auditory and Tactile Masking Noise')
set(gca, 'XScale', 'log'); box off
axis([50 5000 -180 -20])

legend('Auditory Masked')