function [color, newPos] = dfUpdateVisFB(anMsr, rmsMean)
% [color, newPos] = dfUpdateVisFB(anMsr, rmsMean) updates the height and 
% color of the loudness feedback bar that follows each trial in both the 
% SFPerturb and AFPerturb scripts.
%
% The recorded RMS mean from the trial is compared against the baseline
% version that was inputted during dfSetVisFB and the height and color of
% changes based on if the trial recording was within the bounds.

% New Height based on RMS measurement
ratio = (rmsMean - anMsr.vBotAmp)/anMsr.hRMSrange;

newRecHeight = 0.45*ratio;
newDrawHeight = newRecHeight + anMsr.bMar;

% Same X and Y Pos and same X Width...Now with new Height!
newPos = [anMsr.recXSt anMsr.recYSt anMsr.recWidth newRecHeight];

% Compare against the min and max allowable RMS values
if newDrawHeight > anMsr.drawMaxH
    color = 'red';
elseif newDrawHeight < anMsr.drawMinH
    color = 'red';
else
    color = 'green';
end
end