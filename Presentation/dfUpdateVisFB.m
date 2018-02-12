function [color, newPos] = dfUpdateVisFB(anMsr, rmsMean)
%Feedback for the participant on how loud they are. 

%New Height
ratio = (rmsMean - anMsr.vBotAmp)/anMsr.hRMSrange;

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