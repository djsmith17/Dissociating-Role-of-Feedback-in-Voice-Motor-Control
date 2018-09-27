function VisFeedback_testbed
%This function helps to test my visual feedback for RMS present in a
%recording. Currently uses the following functions:

%setPerturbVisualFB.m
%updateVisualFeed.m

curSess     = 'VisFBTest';
recordedRMS = (5+rand(1991,1))/1000;
meanRMS     = mean(recordedRMS);
numTrial    = 4;

targRMS   = 50; %dB just to test
boundsRMS = 3; %+/- dB

close all
[msrStr, annoStr] = dfSetVisFB(curSess, targRMS, boundsRMS);

% pertParadigm(msrStr, annoStr, meanRMS)
% triggerTest(annoStr, numTrial)
endoParadigm(annoStr, numTrial)

close all
end

function pertParadigm(msrStr, annoStr, meanRMS)

%Close the curtains
pause(5); %Let them breathe a sec
set(annoStr.Ready,'Visible','off');

set(annoStr.plus,'Visible','on');
pause(1.0)
set(annoStr.plus,'Visible','off');

set([annoStr.EEE annoStr.visTrig], 'Visible','on');
pause(4.0)
set([annoStr.EEE annoStr.visTrig], 'Visible','off');

[color, newPos] = dfUpdateVisFB(msrStr, meanRMS);

set(annoStr.LoudRec, 'position', newPos);
set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);

set([annoStr.LoudRec, annoStr.fbLines], 'Visible','on'); 
pause(2.0)
set([annoStr.LoudRec, annoStr.fbLines], 'Visible','off');
end

function triggerTest(annoStr, numTrial)
%Close the curtains

pause(1); %Let them breathe a sec
set(annoStr.Ready,'Visible','off');

for i = 1:numTrial
    set(annoStr.plus, 'Visible', 'on');
    pause(1.0)
    set(annoStr.plus, 'Visible', 'off');
    
    set([annoStr.EEE annoStr.visTrig],'Visible','on');
    pause(1.0)
    set([annoStr.EEE annoStr.visTrig],'Visible','off');
    pause(1.0)
end
end

function endoParadigm(annoStr, numTrial)

FBInstr = {'Your Own Voice'; 'Loud White Noise'};

set(annoStr.curSessNote, 'Visible', 'off')
pause(1); %Let them breathe a sec
set(annoStr.Ready,'Visible','off');

for i = 1:numTrial
    curTrial = ['Trial ' num2str(i)];
    AudFB    = FBInstr(mod(i,2)+ 1); % Alternating trials    
    
    set(annoStr.trialT, 'String', curTrial);
    set(annoStr.FBCue, 'String', AudFB);
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'on');
    pause()
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'off');
    
    set(annoStr.plus, 'Visible', 'on');
    pause(1.0)
    set(annoStr.plus, 'Visible', 'off');
    
    set(annoStr.EEE, 'Visible','on');
    pause(1.0)
    set(annoStr.EEE, 'Visible','off');
    pause(1.0)
end

end