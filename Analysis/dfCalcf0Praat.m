function [time, trialsetf0, fsA] = dfCalcf0Praat(dirs, trialset, fs, bTf0b)
%This asks praat to calculate f0 for a given saved wav file. 

resultFolder  = dirs.SavResultsDir;
wavFileLoc    = [resultFolder, '\trialRec.wav'];
txtFileLoc    = [resultFolder, '\pitchCalc.txt'];
[~, numTrial] = size(trialset);

psDir         = dirs.Code;                        %Praat scripting
pbDir         = 'MATLAB-Toolboxes\praatBatching'; %Praat batching

tStep = 0.005; %seconds
fsA   = 1/tStep;

p_fn = fullfile(pbDir, 'praat.exe');
if ~exist(p_fn, 'file')
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pbDir, 'sendpraat.exe');
if ~exist(sp_fn, 'file')
    error('file ''sendpraat.exe'' not found')
end

gt_fn = fullfile(psDir, 'batchcalcf0.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchcalcf0.praat'' not found')
end
 
trialsetf0 = [];
for ii = 1:numTrial
    %Make the audio file in the folder from input data
    audiowrite(wavFileLoc, trialset(:,ii), fs)
    
    curTrial = ii;
    call2 = sprintf('%s praat "execute %s %s %s %f %f', ...
                        sp_fn, ... %sendpraat.exe
                        gt_fn, ... %saved praat script ('generatef0JNDTokens)
                        wavFileLoc, ... %file location of generated wav file
                        txtFileLoc, ...
                        curTrial, ...
                        numTrial ...
                        );

    [s, r] = dos(call2);
    if s ~= 0
        dos([p_fn ' &']);
        [s, r] = dos(call2);
        if s ~= 0
            disp(r)
            error('ERROR: something went wrong')
        end
    end

    praatResult = fopen(txtFileLoc);
    praatScan   = textscan(praatResult, '%f %s');
    fclose(praatResult);
    delete(wavFileLoc);
    delete(txtFileLoc);
    
    [time, trialf0, ~] = PraatPostProcessing(praatScan, bTf0b);
    trialsetf0         = cat(2, trialsetf0, trialf0);
end
end

function [time, f0, meanf0] = PraatPostProcessing(praatScan, bTf0b)
%Expects the table is 2 columns, and the f0 values are in column 2
time   = praatScan{1,1};
f0_str = praatScan{1,2};
recLen = length(f0_str);

praatUndLog = strncmp('--undefined--',f0_str,3);
for ii = 1:recLen
    if praatUndLog(ii) == 1
        f0_str{ii} = num2str(bTf0b);
    end
end
f0   = str2double(f0_str);

meanf0 = mean(f0);
end