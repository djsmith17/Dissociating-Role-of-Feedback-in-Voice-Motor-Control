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
expParam.perCatch      = 1;
expParam.gender        = 'male';
expParam.masking       = 0;
expParam.trialLen      = 4; %Seconds

dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.RecWaveDir    = fullfile(dirs.RecFileDir, 'wavFiles');
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

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

    expParam.cuePause = 1.0;
    expParam.resPause = 2.0;
    
    expParam.targRMS   = 55; %dB %Filler
    expParam.boundsRMS = 3;  %+/- dB
    expParam.win       = 2;  %which monitor? 1 or 2

    %Close the curtains
    [anMsr, H1, H2, fbLines, rec, trigCirc] = dfSetVisFB(expParam.targRMS, expParam.boundsRMS, expParam.win);
    
    DAQin = [];
    for ii = 1:expParam.numTrial
        expParam.curTrial   = ['Trial' num2str(ii)];
        expParam.curSubCond = [expParam.subject expParam.curTrial];
          
        %Level of f0 change based on results from 
        audStimP = []; %dfSetAudapFiles(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType(ii), expParam.trigs(ii,:,1), expParam.stimType);
        
        %Set the OST and PCF functions
        Audapter('ost', expParam.ostFN, 0);
        Audapter('pcf', expParam.pcfFN, 0);
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);   
        
        %Cue to begin trial
        set(H1,'Visible','on');
        pause(expParam.cuePause)
        
        %Phonation Start
        set(H1,'Visible','off');
        set(trigCirc,'Visible','on');
        set(H2,'Visible','on');
        
        fprintf('Running Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');
        Audapter('start');
        
        %This will hold the script for as long as OutputData vector lasts.
        [dataDAQ, time] = s.startForeground;
        
        %Phonation End
        Audapter('stop');
        set(trigCirc,'Visible','off');
        set(H2,'Visible','off'); 

        data = dfSaveRawData(expParam, dirs, p, audStimP, dataDAQ);
        DAQin = cat(3, DAQin, dataDAQ);

        %Grab smooth RMS trace from 'data' structure, compare against baseline
        [color, newPos] = dfUpdateVisFB(anMsr, data.rms(:,1));

        set(rec, 'position', newPos);
        set(rec, 'Color', color); set(rec, 'FaceColor', color);
        set(rec, 'Visible', 'on'); 
        set(fbLines, 'Visible', 'on');  

        pause(expParam.resPause)
        set(fbLines, 'Visible', 'off');
        set(rec, 'Visible', 'off');      
    end
    
    DA.expParam    = expParam;
    DA.dirs        = dirs;
    DA.DAQin       = DAQin;
    DA.audStimP    = audStimP;
    DA.p           = p;

    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject '_DiagAud.mat']);
    save(dirs.RecFileDir, 'DA')
else
    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject '_DiagAud.mat']);
    load(dirs.RecFileDir)
end

niAn = dfAnalysisNIDAQ(DA.expParam, DA.DAQin);

pLimits = [0 4 0 4];
fLimits = [0 4 1 5];
drawDAQsignal(niAn.time, niAn.fSensorC, niAn.fSensorN, niAn.pSensor, niAn.trigs, niAn, pLimits, fLimits, DA.expParam.subject, dirs.SavResultsDir, sv2F)

pLimits = [0 3.5 0 4];
drawDAQcombined(niAn.timeAl, niAn.pSensorAl, niAn.trigs, niAn, pLimits, DA.expParam.subject, dirs.SavResultsDir, sv2F)
end