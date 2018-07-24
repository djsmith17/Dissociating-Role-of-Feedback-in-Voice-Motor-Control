function meanf0 = dfCalcf0PraatSingle(dirs)
%Calculation of pitch using Praat for a single saved wav file.

helperFolder  = dirs.helpers;
tokenDir      = dirs.tokenDir;
wavFileLoc    = dirs.baseTokenFile;
txtFileLoc    = [tokenDir, '\pitchCalc.txt'];
numTrial      = 1;

pbDir         = fullfile(helperFolder, 'praatBatching'); % Praat batching

lwPitchBnd = 75;
upPitchBnd = 300;

p_fn = fullfile(pbDir, 'praat.exe');
if ~exist(p_fn, 'file')
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pbDir, 'sendpraat.exe');
if ~exist(sp_fn, 'file')
    error('file ''sendpraat.exe'' not found')
end

gt_fn = fullfile(pbDir, 'batchcalcf0.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchcalcf0.praat'' not found')
end
 
curTrial = 1;
call2 = sprintf('%s praat "execute %s %s %s %f %f %f %f', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    wavFileLoc, ... %file location for saved wav files
                    txtFileLoc, ...
                    lwPitchBnd, ...
                    upPitchBnd, ... 
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
delete(txtFileLoc);

meanf0      = averagePraatf0(praatScan);
end

function meanf0 = averagePraatf0(praatScan)
%Expects the table is 2 columns, and the f0 values are in column 2

praatf0idx = ~strncmp('--undefined--', praatScan{1,2},3);
praatf0freq = str2double(praatScan{1,2}(praatf0idx));
meanf0 = mean(praatf0freq);
end