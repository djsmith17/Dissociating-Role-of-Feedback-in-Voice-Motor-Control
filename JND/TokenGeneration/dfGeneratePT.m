function PertTokens = dfGeneratePT(dirs, GT, PertFreqs)   
%This expects that you have calculated f0 elsewhere and have already
%determined the spacing in (Hz) for each set of stimuli

subj      = GT.subject;
run       = GT.run;
numTokens = length(PertFreqs);

helperFolder = dirs.helpers;
tokenDir     = dirs.tokenDir;

pbDir     = fullfile(helperFolder, 'praatBatching'); %Praat batching

tokenDir = [tokenDir, '\']; %add a slash to the mic folder
ext = '.wav'; %extension of files

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
 
PertTokens = [];
for ii = 1:numTokens
    targetPert = PertFreqs(ii);
    if GT.xAll(ii) < 0
        pertStr = sprintf('%0.1f', GT.xAll(ii));
    else
        pertStr = sprintf('+%0.1f', GT.xAll(ii));
    end
    targetPertName = [subj run pertStr];

    %Build DOS calls to control praat
    call2 = sprintf('%s praat "execute %s %s %s %s %f %f %f"', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    tokenDir, ... %file location for saved wav files
                    ext, ... %file extension
                    targetPertName, ...
                    targetPert, ... %Frequency to shift to
                    ii, ... %current token being created
                    numTokens ... %Total number of tokens being created
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
    
    [thisToken, ~] = audioread(fullfile(tokenDir, [targetPertName '.wav']));
    PertTokens = cat(1, PertTokens, thisToken');
end
end