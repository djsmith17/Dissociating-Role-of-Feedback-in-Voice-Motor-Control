function [meanJND, lastSetAccu] = dfAnalyzeThresholdJND(UD, varargin)
% modified from Palam
% Received from Ayoub Daliri
HighReversal = max(UD.reversal);
NumTrials = length(UD.response);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'reversals',4)            
            criterion = 'reversals';
            number = varargin{n+1};
            LowReversal = HighReversal - number + 1;
            if LowReversal < 1
                warning('PALAMEDES:invalidOption','You asked for the last %s reversals. There are only %s reversals.',int2str(number),int2str(HighReversal));               
                LowReversal = 1;
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'trials',4)            
            criterion = 'trials';
            number = varargin{n+1};
            LowTrial = NumTrials - number + 1;
            if LowTrial < 1
                warning('PALAMEDES:invalidOption','You asked for the last %s trials. There are only %s trials.',int2str(number),int2str(NumTrials));               
                LowTrial = 1;
            end
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})
        end        
    end            
else
    criterion = 'reversals';
    LowReversal = 3;
end

if strncmpi(criterion,'reversals',4)
    meanJND = sum(UD.xStaircase(UD.reversal >= LowReversal))/(HighReversal-LowReversal+1);
    ind = find(UD.reversal == (LowReversal-1))+1;
    lastSet = UD.catchResponse(ind:end);
    lastSetL = length(lastSet);
    lastSetAccu = 100*sum(lastSet)/lastSetL;
else
    meanJND = sum(UD.xStaircase(LowTrial:NumTrials))/(NumTrials-LowTrial+1);
    lastSetAccu =[];
end