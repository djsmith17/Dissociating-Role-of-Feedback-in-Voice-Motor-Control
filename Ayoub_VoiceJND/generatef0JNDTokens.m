function generatef0JNDTokens(dirs, baseToken)   

tokenDir = dirs.tokenDir;
psDir    = fullfile(dirs.code, presentation);        %Praat scripting
pbDir    = fullfile(dirs.helpers, 'praat batching'); %Praat batching

newmic = [micFolder, '\']; %add a slash to the mic folder
ext = '.wav'; %extension of files

p_fn = fullfile(pbDir, 'praat.exe');
if ~exist(p_fn, 'file')
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pbDir, 'sendpraat.exe');
if ~exist( sp_fn , 'file' )
    error( 'file ''sendpraat.exe'' not found')
end

gp_fn = fullfile(psDir, 'generatef0JNDTokens');
if ~exist(gp_fn, 'file')
    error('file ''generatef0JNDTokens'' not found')
end
    
%Build DOS calls to control praat
call2 = sprintf('%s praat "execute %s %s %s %s %f %f %f %f %f"', ...
                sp_fn, ... %location of sendpraat.exe
                gp_fn, ... %location of praat script wrote (get_fo_EET.praat)
                newmic, ... %where wav files are located aka where you ONE new file is made
                ext, ... %extension of files
                txtmicFolder, ... %destination where .txt files will stay % where you want your output to save (wav file)
                voicingThreshold,  ... %voicing threshold indicated %probably won't need
                pitchFloor, ... %pitch floor indicated %might want to set
                pitchCeiling, ... %pitch ceiling indicated %might want to set
                maxformant, ... %max formant value in % won't need
                numformant ... %number of formants to calculate %won't need
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