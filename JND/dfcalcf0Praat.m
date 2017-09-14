function dfcalcf0Praat(dirs)
%This asks praat to calculate f0 for a given saved wav file. 

tokenDir = dirs.tokenDir;
psDir    = dirs.Code;                        %Praat scripting
pbDir    = 'MATLAB-Toolboxes\praatBatching'; %Praat batching

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

gt_fn = fullfile(psDir, 'batchcalcf0.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchcalcf0.praat'' not found')
end
 
call2 = sprintf('%s praat "execute %s %s %s %s %f"', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    tokenDir, ... %file location for saved wav files
                    ext, ... %file extension
                    targetPertName, ...
                    targetPert ... %Number of Tokens to create
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


PertTokens = [];
for ii = 1:numTokens
    targetPert = PertFreqs(ii);
    targetPertName = ['Cent' num2str(ii)];

    %Build DOS calls to control praat
    call2 = sprintf('%s praat "execute %s %s %s %s %f"', ...
                    sp_fn, ... %sendpraat.exe
                    gt_fn, ... %saved praat script ('generatef0JNDTokens)
                    tokenDir, ... %file location for saved wav files
                    ext, ... %file extension
                    targetPertName, ...
                    targetPert ... %Number of Tokens to create
                    );    


    
    [thisToken, fs] = audioread(fullfile(tokenDir, [targetPertName '.wav']));
    PertTokens = cat(1, PertTokens, thisToken');
end



end