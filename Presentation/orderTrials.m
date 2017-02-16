function trialType = orderTrials(numTrial, per)
%This function organizes trials in sets based on variable 'per'. 
%For example if per = 0.25, then a set of four trials will have 
%three control trials, and one catch trial. The catch trial will always be 
%the second or third trial in the set. This will generate enough sets for
%all trials in 'numTrial'.

r = rem(numTrial,per);

trialType = [];
for ii = 1:(numTrial*per)
    place = round(rand) + 2; %This will break if per!= 0.25 **Fix this**
    set = zeros(1,(1/per));
    set(place) = 1;
    
    trialType = cat(2, trialType, set);
end
end