function [color, newPos] = updateVisualFeed(anMsr, rms)
%Feedback for the participant on how loud they are. 

%This might be a little loose. Confirm if this is ok...
RMS = mean(rms(:,1));

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