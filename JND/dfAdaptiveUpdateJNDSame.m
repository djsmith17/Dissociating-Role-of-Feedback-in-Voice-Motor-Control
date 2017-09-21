function UD = dfAdaptiveUpdateJNDSame(UD, response)
% Modified from Palam
% Received from Ayoub Daliri
trial = length(UD.x);
UD.response(trial) = response;

if trial == 1
    UD.xStaircase(trial) = UD.x(trial);
    if response == 1
        UD.direction = -1;        
    else
        UD.direction = 1;
    end
end
    
if response == 1
    UD.d = UD.d + 1;
    if UD.d == UD.down || max(UD.reversal) < 1
        UD.xStaircase(trial+1) = UD.xStaircase(trial)-UD.stepSizeDown;
        if UD.xStaircase(trial+1) < UD.xMin && strcmp(UD.truncate,'yes')
            UD.xStaircase(trial+1) = UD.xMin + 1; %Cara added the .1 Hz here because if it is truly xMin, the correct answer isn't really different anymore.
        end
        UD.u = 0;
        UD.d = 0;
        UD.reversal(trial) = 0;
        if UD.direction == 1
            UD.reversal(trial) = sum(UD.reversal~=0) + 1;
        else
            UD.reversal(trial) = 0;
        end
        UD.direction = -1;
    else
        UD.xStaircase(trial+1) = UD.xStaircase(trial);
    end    
else
    UD.u = UD.u + 1;
    if UD.u == UD.up || max(UD.reversal) < 1
        UD.xStaircase(trial+1) = UD.xStaircase(trial)+UD.stepSizeUp;
        if UD.xStaircase(trial+1) > UD.xMax && strcmp(UD.truncate,'yes')
            UD.xStaircase(trial+1) = UD.xMax;
        end
        UD.u = 0;
        UD.d = 0;
        UD.reversal(trial) = 0;
        if UD.direction == -1
            UD.reversal(trial) = sum(UD.reversal~=0) + 1;
        else
            UD.reversal(trial) = 0;
        end
        UD.direction = 1;
    else
        UD.xStaircase(trial+1) = UD.xStaircase(trial);
    end    
end    

if strncmpi(UD.stopCriterion,'reversals',4) && sum(UD.reversal~=0) == UD.stopRule
    UD.stop = 1;
end
if strncmpi(UD.stopCriterion,'trials',4) && trial == UD.stopRule
    UD.stop = 1;
end
if ~UD.stop
    UD.x(trial+1) = UD.xStaircase(trial+1);
    if UD.x(trial+1) > UD.xMax
        UD.x(trial+1) = UD.xMax;
    elseif UD.x(trial+1) < UD.xMin
        UD.x(trial+1) = UD.xMin;
    end
    UD.xCurrent = UD.x(trial+1);
end