function trialType = orderTrials(numTrial, per)
%This function organizes trials in sets based on variable 'per'. The first
%trial and the last trial of each set is always a control trial. A random
%trial between the first and last trial will be a catch trial, with the
%remaining trials of the set being control. The length of a set will depend
%on the total number of trials (numTrial), and the percent of catch trials
%(per)in the paradigm. For example if per = 0.25, then a set of trials will
%contain four trials; three control trials, and one catch trial. 
%The catch trial will always be either the second or third trial in the set.  
%This will generate enough sets for all trials in 'numTrial'.

%numTrial: The total number of trials 
%per:      The percent (in decimal) of total trials that are catch trials.

%trialType: A vector of length trialType of 0s and 1s representing the
%order of control and catch trials. 

setSize = 1/per;
options = setSize - 3 ; %Excluding first and last pos in the set and pos 2. 

trialType = [];
numSets = numTrial*per;
if round(numSets) ~= numSets
    disp('ERROR: Give me a whole number of sets')
elseif per == 0
    disp('No trials will be perturbed!')
    trialType = zeros(1, numTrial);
elseif per == 1
    disp('All trials will be perturbed!')
    trialType = ones(1, numTrial);
elseif per == 0.5
    disp('Half of all trials will be pertrubed in a semi-random permutation')
    trialType = zeros(1, numTrial);
    trialType(1:numTrial/2) = 1;
    trialType = trialType(randperm(numTrial));
elseif per < 1.0 && per > 0.5
    disp('ERROR: This only supports percents of 100%, 50%, or lower than 50%')
else
    for ii = 1:numSets
        set = zeros(1, setSize);
        place = round(options*rand) + 2; %Always want to at least start at pos 2    
        set(place) = 1;

        trialType = cat(2, trialType, set);
    end
end  
end