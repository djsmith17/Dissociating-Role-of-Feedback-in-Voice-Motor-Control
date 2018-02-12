function [color, newPos] = dfUpdateVisFB(anMsr, rms)
%Feedback for the participant on how loud they are. 

%Convert rms SPL to dB and take the mean
refSPL  = 0.00002;      %20 micropascals for air
rmsdB   = 20*log10(rms/refSPL);

%There were -Inf in my RMSdB. not exactly sure why, but this fixes the
%trial 1, lack of feeedback
rmsdB(find(rmsdB == -Inf)) = 0;

meanAmp = mean(rmsdB);

%New Height
ratio = (meanAmp - anMsr.vBotAmp)/anMsr.hRMSrange;

newRecHeight = 0.45*ratio;
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