function PertTokens = dfGeneratePT(dirs, GT)   
% dfGeneratePT(dirs, GT) sends praat batching commands to generate 
% individual pitch-shifted tokens from a baseline voice token. Praat shifts
% the f0 of the given baseline recording to a frequency value that is 
% pre-calculated. .wav files are saved with the names of the cent-shift 
% (from baseline) for each each token generated. The token signals are also
% stored in a matrix (output variable 'PertTokens'), which is eventually
% saved as a .mat file (e.g. GT1DRF.mat').

subj          = GT.subject;
run           = GT.run;
allCentShifts = GT.allCentShifts;
pertFreqs     = GT.PertFreqs;
numTokens     = GT.numPertToken;
tokenLenP     = GT.tokenLenP;

helperFolder  = dirs.helpers;
tokenDir      = dirs.TokenDir;
baseFile      = dirs.BaseFile;

pbDir    = fullfile(helperFolder, 'praatBatching'); % Where to find Praat batching

tokenDir = [tokenDir, '\']; % add a slash to the mic folder
ext      = '.wav';          % File extension

p_fn = fullfile(pbDir, 'praat.exe');
if ~exist(p_fn, 'file')
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pbDir, 'sendpraat.exe');
if ~exist(sp_fn, 'file')
    error('file ''sendpraat.exe'' not found')
end

gt_fn = fullfile(pbDir, 'batchf0JNDTokens.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchf0JNDTokens.praat'' not found')
end
 
PertTokens = zeros(numTokens, tokenLenP);
for ii = 1:numTokens
    targetPert = pertFreqs(ii);
    if allCentShifts(ii) < 0
        pertStr = sprintf('%0.1f', allCentShifts(ii));
    else
        pertStr = sprintf('+%0.1f', allCentShifts(ii));
    end
    targetPertName = [subj run pertStr];

    %Build DOS calls to control praat
    call2 = sprintf('%s praat "execute %s %s %s %s %f %f %f"', ...
                    sp_fn, ...      % sendpraat.exe
                    gt_fn, ...      % batchf0JNDTokens.praat
                    tokenDir, ...   % file location for saved wav files
                    ext, ...        % file extension
                    targetPertName, ...
                    targetPert, ... % Target Frequency Shift
                    ii, ...         % Current token being created
                    numTokens ...   % Total number of tokens
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
    
    [thisToken, ~]    = audioread(fullfile(tokenDir, [targetPertName '.wav']));
    PertTokens(ii, :) = thisToken;
end
end