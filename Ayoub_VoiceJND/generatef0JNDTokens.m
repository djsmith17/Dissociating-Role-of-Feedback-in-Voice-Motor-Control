function tokens = generatef0JNDTokens(dirs, baseToken)   

tokenDir = dirs.tokenDir;
psDir    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Analysis';        %Praat scripting
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

gt_fn = fullfile(psDir, 'batchf0JNDTokens.praat');
if ~exist(gt_fn, 'file')
    error('file ''batchf0JNDTokens.praat'' not found')
end
    
%Build DOS calls to control praat
call2 = sprintf('%s praat "execute %s %f %f"', ...
                sp_fn, ... %sendpraat.exe
                gt_fn, ... %saved praat script ('generatef0JNDTokens)
                tokenDir, ... %file location for saved wav files
                ext ... %file extension
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
end