function meanf0 = dfcalcf0Praat(dirs)
%This asks praat to calculate f0 for a given saved wav file. 

tokenDir = dirs.tokenDir;
psDir    = dirs.JNDTG;                       %Praat scripting
pbDir    = 'MATLAB-Toolboxes\praatBatching'; %Praat batching

tokenDir   = [tokenDir, '\']; % add a slash to the mic folder
ext        = '.wav';          % extension of files
txtFileLoc = [tokenDir, '\pitchCalc.txt'];

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
 
call2 = sprintf('%s praat "execute %s %s %s %s', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    tokenDir, ... %file location for saved wav files
                    ext, ... %file extension
                    txtFileLoc ...
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
meanf0      = averagePraatf0(praatScan);

fclose(praatResult);
delete(txtFileLoc);
end

function meanf0 = averagePraatf0(praatScan)
%Expects the table is 2 columns, and the f0 values are in column 2

praatf0idx = ~strncmp('--undefined--', praatScan{1,2},3);
praatf0freq = str2double(praatScan{1,2}(praatf0idx));
meanf0 = mean(praatf0freq);
end