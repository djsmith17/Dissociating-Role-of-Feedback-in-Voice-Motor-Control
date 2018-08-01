function dfShowPraatSpect(dirs, curSess, trial)
%This asks praat to calculate f0 for a given saved wav file. 

helperFolder   = dirs.helpers;
wavFolder      = dirs.RecWaveDir;
wavFileLoc     = fullfile(wavFolder, [curSess trial '_micIn.wav']);

if ~exist(wavFileLoc, 'File')
    fprintf('Wav file does not exist')
    return
end

pbDir         = fullfile(helperFolder, 'praatBatching'); % Praat batching

p_fn = fullfile(pbDir, 'praat.exe');
if ~exist(p_fn, 'file')
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pbDir, 'sendpraat.exe');
if ~exist(sp_fn, 'file')
    error('file ''sendpraat.exe'' not found')
end

gt_fn = fullfile(pbDir, 'batchShowPraatSpect.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchShowPraatSpect.praat'' not found')
end
 
call2 = sprintf('%s praat "execute %s %s %s %f %f', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    wavFileLoc ... %file location of generated wav file
                    );

[s, ~] = dos(call2);
if s ~= 0
    dos([p_fn ' &']);
    [s, r] = dos(call2);
    if s ~= 0
        disp(r)
        error('ERROR: something went wrong')
    end
end
end