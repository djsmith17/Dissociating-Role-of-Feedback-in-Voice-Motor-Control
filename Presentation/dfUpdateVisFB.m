function [color, newRecMeas, loudResult] = dfUpdateVisFB(msrStr, sampRMS)
% [color, newPos, loudResult] = dfUpdateVisFB(msrStr, sampRMS) updates the 
% height and color of the loudness bar that is present at the end of trials 
% as visual feedback to the participant.
% This script is updating the annotation elements that are created by 
% the function 'presentation\dfSetVisFB.m'
%
% The mean RMS from a sample/trial is compared against the baseline
% sample that was inputted during dfSetVisFB and the height and color of
% changes based on if the trial recording was within the bounds.
%
% INPUTS
% msrStr:  Struc of annotation measurements. The collection of the 
%          positions for all the annotations on the figure. Annotations are
%          just objects on the figure that can be moved around. This struct
%          is generated from 'presentation\dfSetVisFB.m'
% sampRMS: Mean RMS of sample trial recording
%
% OUTPUTS
% color:      String ('Red', 'Green'). Color update to the loudness bar
% newRecMeas: Vector [XPos YPos Width Height] of the updated position for the 
%             loudness bar
% loudResult: Int value for the outcome of the loudness comparison (0 = Just 
%             Right, -1 = Too Low, 1 = Too High)
%
% Author: Dante J Smith
% Last Update: 10/28/2021

% New Height based on RMS measurement
ratio = (sampRMS - msrStr.vBotAmp)/msrStr.hRMSrange;

newRecHeight = 0.45*ratio;
newDrawHeight = newRecHeight + msrStr.bMar;

% Same X/Y pos and same width...now with new height!
newRecMeas = [msrStr.recXSt msrStr.recYSt msrStr.recWidth newRecHeight];

% Compare against the min/max allowable RMS values
if newDrawHeight > msrStr.drawMaxH
    color = 'red';
    loudResult = 1;  # Too loud
elseif newDrawHeight < msrStr.drawMinH
    color = 'red';
    loudResult = -1; # Too soft
else
    color = 'green';
    loudResult = 0;  # Just right
end
end
