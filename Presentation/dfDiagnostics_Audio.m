function dfDiagnostics_Audio(varargin)
%A quick test of the force sensors before running the actual experiment.
%This makes sure that the sensors are working they should be and we can
%continue with the experiment. Eventually this will also include the
%pressure sensor. 

%This script calls the following (4) functions:
%dfDirs.m
%initNIDAQ.m
%dfMakePertSignal.m
%drawDAQsignal.m
close all;

if isempty(varargin)
    numTrial = 10; 
else
    numTrial = varargin{1};
end

collectNewData         = 1; %Boolean
sv2F                   = 1; %Boolean

expParam.project       = 'Diagnostics_Audio';
expParam.expType       = 'Auditory Perturbation_Perceptual';
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.numTrial      = numTrial; %Experimental trials = 40
expParam.trialLen      = 4; %Seconds
expParam.perCatch      = 1;
expParam.gender        = 'male';

dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject '_DiagAud.mat']);

if collectNewData == 1
    %Paradigm Configurations
    expParam.sRate              = 48000;
    expParam.downFact           = 3;
    expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.frameLen           = 96;  % Before downsampling
    expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'
    
    %Set up Audapter
    Audapter('deviceName', expParam.audioInterfaceName);
    Audapter('setParam', 'downFact', expParam.downFact, 0);
    Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
    Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
    p = getAudapterDefaultParams(expParam.gender);

    [s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev3');
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;
    
    %Set up OST and PCF Files
    expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
    expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

    expParam.resPause = 1;

    DAQin = [];
    for ii = 1:expParam.numTrial
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);        
        
        fprintf('Running Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');
        Audapter('start');
        
        %This will hold the script for as long as OutputData vector lasts.
        [data_DAQ, time] = s.startForeground;
        
        %Phonation End
        Audapter('stop');

        DAQin = cat(3, DAQin, data_DAQ);

        pause(expParam.resPause)      
    end
    
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = DAQin;

    save(dirs.RecFileDir, 'NSD')
else
    load(dirs.RecFileDir)
end

niAn = dfAnalysisNIDAQ(NSD.expParam, NSD.DAQin);

pLimits = [0 4 0 4];
fLimits = [0 4 1 5];
drawDAQsignal(niAn.time, niAn.fSensorC, niAn.fSensorN, niAn.pSensor, niAn.trigs, niAn, pLimits, fLimits, NSD.expParam.subject, dirs.SavResultsDir, sv2F)

pLimits = [0 3.5 0 4];
drawDAQcombined(niAn.timeAl, niAn.pSensorAl, niAn.trigs, niAn, pLimits, NSD.expParam.subject, dirs.SavResultsDir, sv2F)
end