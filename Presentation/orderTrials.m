function trialType = orderTrials(numTrial, per)
%This function organizes trials in sets based on variable 'per'. 
%For example if per = 0.25, then a set of four trials will have 
%three control trials, and one catch trial. The catch trial will always be 
%the second or third trial in the set. This will generate enough sets for
%all trials in 'numTrial'.

setSize = 1/per;
options = setSize -2 ; %Excluding first and last position in the set. 

trialType = [];
numSets = numTrial*per;
if round(numSets) ~= numSets
    disp('ERROR: Give me a whole number of sets')
elseif per >= 0.5
    disp('ERROR: This pattern does not support 50% catch trials')
else
    for ii = 1:numSets
        set = zeros(1, setSize);
        place = round(rand) + options ;    
        set(place) = 1;

        trialType = cat(2, trialType, set);
    end
end  
end