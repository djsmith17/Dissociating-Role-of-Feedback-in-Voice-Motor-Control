   
newmic = [micFolder, '\']; %add a slash to the mic folder
ext = '.wav'; %extension of files


%% make sure you have the scripts you need for praat interface to work
pi_dir = 'praat_interface' ;
if ~exist(pi_dir , 'dir')
    error('folder ''praat_interface'' not found in directory')
end

p_fn = fullfile(pi_dir , 'praat.exe') ;
if ~exist(p_fn , 'file' )
    error('file ''praat.exe'' not found')
end

sp_fn = fullfile(pi_dir , 'sendpraat.exe') ;
if ~exist( sp_fn , 'file' )
    error( 'file ''sendpraat.exe'' not found' )
end
gp_fn = fullfile( pi_dir , 'get_f0_EET.praat' ) ;
if ~exist( gp_fn , 'file' )
    error( 'file ''get_f0_EET.praat'' not found' )
end
    
%% uncomment to run praat script
%Build DOS calls to control praat
call2 = sprintf( '%s praat "execute %s %s %s %s %f %f %f %f %f"' , ...
    fullfile( pwd , sp_fn ) , ... %location of sendpraat.exe
    fullfile( pwd , gp_fn ) , ... %location of praat script wrote (get_fo_EET.praat)
    newmic , ... %where wav files are located aka where you ONE new file is made
    ext , ... %extension of files
    txtmicFolder, ... %destination where .txt files will stay % where you want your output to save (wav file)
    voicingThreshold,  ... %voicing threshold indicated %probably won't need
    pitchFloor, ... %pitch floor indicated %might want to set
    pitchCeiling, ... %pitch ceiling indicated %might want to set
    maxformant, ... %max formant value in % won't need
    numformant ... %number of formants to calculate %won't need
    ) ;
    
    
    
[s, r] = dos(call2);
if s ~= 0
    dos([fullfile(pwd, p_fn) ' &']);
    [s, r] = dos(call2);
    if s ~= 0
        disp(r)
        error('ERROR: something went wrong')
    end
end