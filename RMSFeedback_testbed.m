function RMSFeedback_testbed
%This function helps to test my visual feedback for RMS present in a
%recording. Currently uses the following functions:

%setPerturbVisualFB.m
%updateVisualFeed.m

recordedRMS = (5+rand(1991,1))/1000;

targRMS   = 50; %dB just to test
boundsRMS = 3; %+/- dB
win = 2;

close all
[anMsr, H1, H2, fbLines, rec, trigCirc] = dfSetVisFB(targRMS, boundsRMS, win);

pause(1.0)
set(H1,'Visible','off');

[color, newPos] = dfUpdateVisFB(anMsr, recordedRMS);

set(rec, 'position', newPos);
set(rec, 'Color', color); set(rec, 'FaceColor', color);
set(rec, 'Visible','on'); 
set(fbLines, 'Visible', 'on'); 
set(trigCirc, 'Visible', 'on')

pause(2.0)
set(trigCirc, 'Visible', 'off')
pause(2.0)
set(trigCirc, 'Visible', 'on')
pause(2.0)
set(fbLines, 'Visible','off');
set(rec, 'Visible','off');
set(trigCirc, 'Visible', 'off')
close all
end