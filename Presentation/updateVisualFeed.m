function [color, newPos] = updateVisualFeed(anMsr, rms)
%Feedback for the participant on how loud they are. 

%Convert rms SPL to dB and take the mean
refSPL  = 0.00002;      %20 micropascals for air
rmsdB   = 20*log10(rms/refSPL);
rmsMean = mean(rmsdB);

%New Height
newRecHeight = 0.55; %This needs to be a calculation of current RMS against maximal RMS (anMsr)
newDrawHeight = newRecHeight + anMsr.bMar;

newPos = [anMsr.recXSt anMsr.recYSt anMsr.recWidth newRecHeight];

if newDrawHeight > anMsr.drawMaxH
    color = 'red';
elseif newDrawHeight < anMsr.drawMinH
    color = 'red';
else
    color = 'green';
end
end