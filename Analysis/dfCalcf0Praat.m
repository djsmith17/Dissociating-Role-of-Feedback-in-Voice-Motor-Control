function trialsetf0 = dfCalcf0Praat(dirs, trialset, fs)
%This asks praat to calculate f0 for a given saved wav file. 

resultFolder = dirs.SavResultsDir;
wavFileLoc   = [resultFolder, '\trialRec.wav'];
txtFileLoc   = [resultFolder, '\pitchCalc.txt'];
psDir        = dirs.Code;                        %Praat scripting
pbDir        = 'MATLAB-Toolboxes\praatBatching'; %Praat batching
numTrial     = length(trialset);

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
    trialf0     = PraatPostProcessing(praatScan);
    trialsetf0  = cat(2, trialsetf0, trialf0);

    fclose(praatResult);
    delete(wavFileLoc);
    delete(txtFileLoc);
end
end

function meanf0 = PraatPostProcessing(praatScan)
%Expects the table is 2 columns, and the f0 values are in column 2

praatf0idx = ~strncmp('--undefined--', praatScan{1,2},3);
praatf0freq = str2double(praatScan{1,2}(praatf0idx));
meanf0 = mean(praatf0freq);
end