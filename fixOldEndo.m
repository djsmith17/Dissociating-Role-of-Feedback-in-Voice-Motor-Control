dirs = 'C:\Users\djsmith\Desktop\SFL1';

subject = 'DRF_ENP1';
run     = 'SFL1';
trials  = 1:5;

allDRF      = [];
trialType   = [];
sigs        = [];
trigs       = [];
for i = trials
    load(fullfile(dirs, [subject run 'Trial' num2str(i) 'DRF.mat']));

    allDRF = cat(1, allDRF, DRF);
    curExpParam = DRF.expParam;
    
    trialType = cat(1, trialType, curExpParam.trialType);
    sigs      = cat(2, sigs, curExpParam.sigs);
    trigs     = cat(1, trigs, curExpParam.trigs);
end

numTrial = length(trials);
curSess = [subject run];

for j = trials
   curDRF      = allDRF(j);
   curExpParam = curDRF.expParam;
   
   curExpParam.run      = run;
   curExpParam.numTrial = numTrial;
   curExpParam.curTrial = ['Trial' num2str(j)];
   curExpParam.curSess  = curSess;
   curExpParam.curSessTrial = [curSess curExpParam.curTrial];
    
    
   curExpParam.trialType = trialType;
   curExpParam.sigs      = sigs;
   curExpParam.trigs     = trigs;
   
   curDRF.expParam = curExpParam;
   
   DRF = curDRF;
   save(fullfile(dirs, [subject run 'Trial' num2str(j) 'DRF.mat']), 'DRF');
end