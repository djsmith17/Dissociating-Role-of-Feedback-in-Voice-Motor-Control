function RMSFeedback_testbed
%This function helps to test my visual feedback for RMS present in a
%recording. Currently uses the following functions:

%setPerturbVisualFB.m
%updateVisualFeed.m

recordedRMS = (5+rand(1991,1))/1000;
numTrial = 1;

targRMS   = 50; %dB just to test
boundsRMS = 3; %+/- dB

close all
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(targRMS, boundsRMS);

% pertParadigm(anMsr, H1, H2, H3, fbLines, rec, trigCirc, recordedRMS)
triggerTest(anMsr, H1, H2, H3, fbLines, rec, trigCirc, numTrial)

close all
end

function pertParadigm(anMsr, H1, H2, H3, fbLines, rec, trigCirc, recordedRMS)

%Close the curtains
pause(5); %Let them breathe a sec
set(H3,'Visible','off');

set(H1,'Visible','on');
pause(1.0)
set(H1,'Visible','off');

set(H2,'Visible','on');
set(trigCirc,'Visible','on');
pause(4.0)
set(H2,'Visible','off');
set(trigCirc,'Visible','off');

[color, newPos] = dfUpdateVisFB(anMsr, recordedRMS);

set(rec, 'position', newPos);
set(rec, 'Color', color); set(rec, 'FaceColor', color);

set(rec, 'Visible','on'); 
set(fbLines, 'Visible', 'on'); 
pause(2.0)
set(rec, 'Visible','off');
set(fbLines, 'Visible','off');
end

function triggerTest(anMsr, H1, H2, H3, fbLines, rec, trigCirc, numTrial)
%Close the curtains
pause(1); %Let them breathe a sec
set(H3,'Visible','off');

for i = 1:numTrial
    set([H2 trigCirc],'Visible','on');
    pause(1.0)
    set([H2 trigCirc],'Visible','off');
    pause(1.0)
end
end