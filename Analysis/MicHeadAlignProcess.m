classdef MicHeadAlignProcess
    % MicHeadAlignProcess(analysisVar, trialVar) is a class which handles
    % the time-alignment between audio signals collected in individual 
    % trials of sensorimotor experiments. 
    % 
    % This calculates delays between the Audapter and NIDAQ microphone
    % recordings, as well as the delay between the Microphone and
    % Headphones in the Audapter recording
    %
    % analysisVar is a structure of variables from the overall experiment
    % and expects the following (3) fields.
    % -expType   : The Experiment Name (e.g. 'Auditory Perturbation_Perceptual')
    % -sRate     : Sampling rate of recordings
    % -frameLen  : Frame length of recordings
    %
    % trialVar is a structure of variables from the specific trial being
    % processed and expects the following (14) fields.
    % -AudFB
    % -rawMic
    % -rawHead
    % -rms
    % -auTrigs
    % -fsNi
    % -trialTimeNI
    % -expTrigsNI
    % -timeNI
    % -rawMicNI
    % -pressureNI
    % -pressureTrigs
    % -presLagTimes
    % -presRiseTimes
    
    properties
        expType     % Experiment Type
        AudFB       % Auditory Feedback provided
        rawMic      % Raw microphone signal
        rawHead     % Raw headphone signal
        rms         % rms from microphone signal
        fs          % sampling rate of audio recordings
        frameLen    % Frame rate of recording (After downsampling)
        trialLen    % length (points) of raw microphone/headphone recording
        trialTime   % length (time) of raw microphone/headphone recording
        time        % time vector of the recording
        auTrigs     % trigger points from Audapter recording
        auTrigsAuNi % trigger points from Audapter recording aligned with NIDAQ triggers (Only used for quality check)
        auTimesAuNi % time points from Aduapter recording aligned with NIDAQ triggers (Only used for quality check)
        prePertVoicingTime % amount of vocalization pre-perturbation
        rawMicDS    % Raw microphone signal, downsampled
        
        frameDel    % Expected frame-delay between Microphone and Headphones (Audapter)
        rmsThresh
        
        fsNI
        trialTimeNI
        expTrigsNI
        timeNI
        rawMicNI
        pressureNI
        pressureTrigs
        presLagTimes
        presRiseTimes
        pressureEndActions
        
        numSamp
        
        rmsVoiceFrame
        voiceOnsetInd
        voiceOnsetT
        
        preVOnsetRMS
        
        AuNIDelay    % Delay between the Audapter and NIDAQ recordings in time
        AuNIDelayP   % Delay between the Audapter and NIDAQ recordings in points (Au)
        adjustedDelay
        
        AuMHDelay    % Delay between the Microphone and Headphone recordings in time
        AuMHDelayP   % Delay between the Microphone and Headphone recordings in points (Au)
        
        analysisSec
        analysisPoints
        analysisFrames
        
        voiceOnsetLate
        breakOccured
        tooShort
        
        processedMic
        processedHead
        
        saveT
        saveTmsg
    end
    
    methods
        function obj = MicHeadAlignProcess(analysisVar, trialVar)
            % obj = MicHeadAlignProcess(analysisVar, trialVar) aligns 
            
            % Experimental Variables
            obj.expType   = analysisVar.expType;
            obj.AudFB     = trialVar.AudFB;
            
            % Audapter Recorded Signals and Variables
            obj.rawMic    = double(trialVar.rawMic);
            obj.rawHead   = double(trialVar.rawHead);
            obj.rms       = double(trialVar.rms);
            obj.fs        = analysisVar.sRate;
            obj.frameLen  = analysisVar.frameLen;
            obj.trialLen  = length(obj.rawMic);
            obj.trialTime = obj.trialLen/obj.fs;
            obj.time      = linspace(0, obj.trialTime, obj.trialLen);
            obj.auTrigs   = trialVar.auTrigs;
            obj.frameDel  = 7;
            obj.rmsThresh = 0.011;
            
            % NIDAQ Recorded Signals and Variables
            obj.fsNI          = trialVar.fsNI;
            obj.trialTimeNI   = trialVar.trialTimeNI;
            obj.expTrigsNI    = trialVar.expTrigsNI;
            obj.timeNI        = trialVar.timeNI;
            obj.rawMicNI      = trialVar.rawMicNI;
            obj.pressureNI    = trialVar.pressureNI;
            obj.pressureTrigs = trialVar.pressureTrigs;
            obj.presLagTimes  = trialVar.presLagTimes;
            obj.presRiseTimes = trialVar.presRiseTimes;
            obj.pressureEndActions = obj.pressureTrigs + obj.presRiseTimes;
            
            % Define the output NumSamp
            obj.numSamp = obj.trialTimeNI*obj.fs;

            % Identify Voice Onset in the Audapter Microphone (rms) signal
            obj.rmsVoiceFrame = find(obj.rms > obj.rmsThresh);
            obj.voiceOnsetInd = (obj.rmsVoiceFrame(1) - obj.frameDel)*obj.frameLen;
            if obj.voiceOnsetInd <= 0 % Can't have index <= 0
                obj.voiceOnsetInd = 1;
            end
            obj.voiceOnsetT   = obj.time(obj.voiceOnsetInd);
            
            % Evaluate the pre-voice onset rms (used to identify voice breaks)
            obj.preVOnsetRMS = evalPreVoiceRMS(obj);
            
            % Find the delay between NIDAQ recording and Audapter recording
            obj.rawMicDS   = resample(obj.rawMic, obj.fsNI, obj.fs);         % Downsample the Audapter recording
            obj.AuNIDelay  = obj.xCorrTimeLag(obj.rawMicDS, obj.rawMicNI, obj.fsNI); % Perform xCorr between NIDAQ and Audapter. Expect that NIDAQ leads Audapter
            obj.AuNIDelay  = round(obj.AuNIDelay, 3); % Round to the nearest 1ms
            obj.AuNIDelayP = obj.AuNIDelay*obj.fs;    % Convert to points

            % Adjust Triggers against NIDAQ if a Laryngeal Pert Exp.
            % Otherwise adjust based on VoiceOnset, which happens during the PSR
            if strcmp(obj.expType(1:3), 'Som')
                obj.adjustedDelay = obj.AuNIDelayP + round(obj.presLagTimes(1)*obj.fs);
            else
                obj.adjustedDelay = obj.voiceOnsetInd;
            end
            obj.auTrigsAuNi = obj.auTrigs + obj.adjustedDelay;
            obj.auTimesAuNi = obj.time(obj.auTrigsAuNi);
            
            % Amount of Time Vocalizing from Voice Onset to perturbation
            % trigger
            obj.prePertVoicingTime = obj.auTimesAuNi(1) - obj.voiceOnsetT;

            % Aim to section audio at 0.5s pre-onset to 1.0s post-offset.
            preOn   = 0.5*obj.fs;
            postOff = 1.0*obj.fs;

            % Audio points on either side of the perturbation period.
            obj.analysisSec(1) = obj.auTrigsAuNi(1) - preOn;   % Where to start the Analysis period
            obj.analysisSec(2) = obj.auTrigsAuNi(2) + postOff; % Where to end the Analysis period
            if obj.analysisSec(2) > obj.trialLen
                obj.analysisSec(2) = obj.trialLen;
            end    
            obj.analysisPoints = obj.analysisSec(1):obj.analysisSec(2);
            obj.analysisFrames = round(obj.analysisSec(1)/obj.frameLen):round(obj.analysisSec(2)/obj.frameLen);

            % Check the voice onset time against when we want to start analyzing data
            obj.voiceOnsetLate = obj.analysisSec(1) < obj.voiceOnsetInd;

            % Identify if there were any Voice Breaks
            [obj.breakOccured, breakMsg] = obj.identifyVoiceBreak;

            % Find the delay between Audapter Headphone and Microphone
            if obj.AudFB == 2 % No Headphone Out
                obj.AuMHDelay = (obj.frameLen*(obj.frameDel-1))/obj.fs;
            else
                prePertPer = obj.analysisSec(1): obj.auTrigsAuNi(1); % 500ms preperturbation
                obj.AuMHDelay = obj.xCorrTimeLag(obj.rawHead(prePertPer), obj.rawMic(prePertPer), obj.fs);   % Expect Mic leads Head
            end
            obj.AuMHDelay  = round(obj.AuMHDelay, 3); %Round to nearest 1ms
            obj.AuMHDelayP = obj.AuMHDelay*obj.fs; % Convert to points

            %%%%%ADJUSTING LENGTHS OF MIC/HEAD BASED ON DELAYS
            % Align the Microphone and Headphones
            if obj.AuMHDelayP >= 0
                micAuAl  = obj.rawMic(1:(end-obj.AuMHDelayP));
                headAuAl = obj.rawHead((obj.AuMHDelayP+1):end);
            else
                micAuAl  = obj.rawMic;
                headAuAl = obj.rawHead;
            end

            % Adjust for delay between Audapter and NIDAQ
            if obj.adjustedDelay > 0 % As long as the delay is non 0
                micAuNi  = micAuAl(obj.adjustedDelay:end);
                headAuNi = headAuAl(obj.adjustedDelay:end);
            else
                micAuNi  = micAuAl;
                headAuNi = headAuAl;
            end

            % If for whatever reason the audio signal is too short, zero-pad it
            obj.tooShort = 0;
            if length(micAuNi) < obj.numSamp
                diffLen = obj.numSamp - length(micAuNi);
                micAuNi  = [micAuNi; zeros(diffLen, 1)];
                headAuNi = [headAuNi; zeros(diffLen, 1)];
                obj.tooShort = 1;
            end

            % Apply logic for saving a trial
            % Include message for reason saved/not saved
            if obj.voiceOnsetLate
                saveT    = 0;  
                saveTmsg = 'Participant started too late!!';
            elseif obj.breakOccured
                saveT    = 0;
                saveTmsg = 'Participant had a voice break!!';
            elseif obj.tooShort
                if strcmp(obj.expType(1:3), 'Aud')
                    saveT    = 1;
                    saveTmsg = 'Everything is ok, but recording should be longer';
                else
                    saveT    = 0;
                    saveTmsg = 'Recording not long enough';
                end
            else
                saveT    = 1;
                saveTmsg = 'Everything is good';
            end

            % Grab the full numSamp so they can be concatenated cleanly
            obj.processedMic  = micAuNi(1:obj.numSamp);
            obj.processedHead = headAuNi(1:obj.numSamp);

            obj.saveT    = saveT;    % Save trial or no?
            obj.saveTmsg = saveTmsg; % Reason, if any the trial was thrown out
        end
        
        function threshIdx = calcEnvelope(obj, audio, fs)
            % Method for identifying voice onset without the knowing what
            % the pre-voice rms in the signal is thresholded at
            
            envThresh = 0.30; 
            cutoffF   = 40;
            [B, A] = butter(4, cutoffF/(fs/2));

            % Envelope the signal by low-pass filtering (change in amplitude/time ~RMS)
            env = filter(B, A, abs(audio)); 
            
            % Largest peak in the envelope theoretically occurs during voicing
            maxPeak = max(env);

            % Find values that are within threshold of max 'voicing' value
            threshIdx = find(env > envThresh*maxPeak);
        end
        
        function preVOnsetRMS = evalPreVoiceRMS(obj)

            preVOnsetTime   = 0.05;            % 50ms before voice onset
            VOnsetFrame     = floor(obj.voiceOnsetInd/obj.frameLen);
            preVOnsetFrames = floor(preVOnsetTime*obj.fs/obj.frameLen);

            preVoiceRange = (-preVOnsetFrames:0)+VOnsetFrame;
            if sum(preVoiceRange <= 0) > 0
                preVOnsetRMS = obj.rmsThresh;
            else
                preVOnsetRMS = obj.rmsThresh; %mean(rms(preVoiceRange));
            end
        end
        
        function timeLag = xCorrTimeLag(obj, sig1, sig2, fs)
            % xCorrTimeLag(sig1, sig2, fs) calculates the lag between two (seemingly) 
            % identical time based signals. 
            %
            % if timeLag is negative, then sig1 leads sig2. 
            % if timeLag is positive, then sig1 lags sig2.

            % Simple crosscorrelation between two signals
            % Finds the largest peak of the result
            [r, lags]    = xcorr(sig1, sig2);
            r(lags<0) = 0;
            [~, peakInd] = max(r);
            maxLag       = lags(peakInd);
            timeLag      = maxLag/fs;
        end
        
        function [timeSet, delaySet] = MHdelayChunked(obj, sig1, sig2, fs)

            numS = length(sig1);
            chunkL = 0.05;
            chunkP = fs*chunkL;
            numChunk = floor(numS/chunkP);

            timeSet  = zeros(numChunk, 1);
            delaySet = zeros(numChunk, 1);
            for ii = 1:numChunk
                set = (1:chunkP) + (ii-1)*chunkP;

                timeChunk = set(1)/fs;
                sig1Chunk = sig1(set);
                sig2Chunk = sig2(set);

                delay = xCorrTimeLag(sig1Chunk, sig2Chunk, fs);
                timeSet(ii)  = timeChunk;
                delaySet(ii) = delay*1000;
            end
        end
        
        function [breakOccured, breakMsg] = identifyVoiceBreak(obj)

            postVOnsetT     = 0.5; %500ms post VO
            postVOnsetFrame = postVOnsetT*obj.fs/obj.frameLen;
            voiceOnsetFrame = ceil(obj.voiceOnsetInd/obj.frameLen);
            postVOnsetFrames = (0:postVOnsetFrame) + voiceOnsetFrame;

            analysisPerFO  = obj.rms(obj.analysisFrames) < obj.preVOnsetRMS; %Analysis Frame
            postVOnsetFO   = obj.rms(postVOnsetFrames) < obj.preVOnsetRMS;   %Post Voice Onset

            breakTol             = 0.1; % Voice Break Tolerance; 100ms
            breakTolFrame        = breakTol*obj.fs/obj.frameLen;
            breakOccuredAnalysis = sum(analysisPerFO) > breakTolFrame; % Last longer than break tolerance
            breakOccuredPostVO   = sum(postVOnsetFO) > breakTolFrame; % Last longer than break tolerance

            if breakOccuredAnalysis
                breakOccured = 1;
                breakMsg     = 'Break During Analysis';
            elseif breakOccuredPostVO
                breakOccured = 1;
                breakMsg     = 'Break Following VO, Possibly MisIdentified VO';
            else
                breakOccured = 0;
                breakMsg     = '';
            end
        end
        
        function drawPreProcessDiagnostic(obj)

        if strcmp(obj.expType, 'Auditory Perturbation_Perceptual')
            pertStr = 'f0 Shift (cents)';
        else
            pertStr = 'Pressure (psi)';
        end
            
        plotPos = [-1280 200];
        plotDim = [1200 900];
        lineThick = 2.5;
        voiceOnsetColor = 'm';

        % Time Bounds
        auTimeRange = [0 6];
        niTimeRange = auTimeRange - obj.AuNIDelay;

        % Inflation//Deflation Properties
        PertColor = [0 0 0];
        PresColor = [255, 66, 0]/255; %Sunburnt Orange;
        InfColor = 'g';
        DefColor = 'r';
        
        % Setup the Figure parameters
        MHFig = figure('Color', [1 1 1]);
        set(MHFig, 'Position', [plotPos plotDim],'PaperPositionMode','auto')

        ha = tight_subplot(3,1,[0.11 0.05],[0.08 0.08],[0.05 0.05]);

        % Raw NIDAQ Microphone
        axes(ha(1))
        plot(obj.timeNI, obj.rawMicNI)
        hold on
        plot([obj.expTrigsNI(1) obj.expTrigsNI(1)], [-2 2], 'Color', PertColor, 'LineStyle', '-', 'LineWidth', lineThick)
        hold on
        plot([obj.expTrigsNI(2) obj.expTrigsNI(2)], [-2 2], 'Color', PertColor, 'LineStyle', '-', 'LineWidth', lineThick)
        axis([niTimeRange min(obj.rawMicNI) max(obj.rawMicNI)])
        title('Raw NIDAQ Microphone')
        
        % Converted Pressure Recording
        yyaxis right
        plot(obj.timeNI, obj.pressureNI, 'Color', PresColor, 'LineWidth', lineThick)
        
        hold on
        plot([obj.pressureTrigs(1) obj.pressureTrigs(1)], [-600 600], 'Color', InfColor, 'LineStyle', '-', 'LineWidth', lineThick)
        hold on
        plot([obj.pressureEndActions(1) obj.pressureEndActions(1)], [-600 600], 'Color', InfColor, 'LineStyle', '--', 'LineWidth', lineThick)
        
        hold on
        plot([obj.pressureTrigs(2) obj.pressureTrigs(2)], [-600 600], 'Color', DefColor, 'LineStyle', '-', 'LineWidth', lineThick)
        hold on
        plot([obj.pressureEndActions(2) obj.pressureEndActions(2)], [-600 600], 'Color', DefColor, 'LineStyle', '--', 'LineWidth', lineThick)       
        
        ylabel(pertStr, 'Color', PresColor)
        axis([niTimeRange -0.2 5.5])
        box off  

        xlabel('NIDAQ Time (s)')
        set(gca,'FontName', 'Arial',...
                'FontSize', 14,...
                'FontWeight','bold',...
                'YColor', PresColor)

        % Raw Audapter Microphone
        axes(ha(2))
        plot(obj.time, obj.rawMic)
        hold on
        plot([obj.voiceOnsetT obj.voiceOnsetT], [-2 2], 'Color', voiceOnsetColor, 'LineStyle', '--', 'LineWidth', lineThick)
        hold on
        plot([obj.time(obj.auTrigs(1)) obj.time(obj.auTrigs(1))], [-2 2], 'k--', 'LineWidth', lineThick)
        hold on
        plot([obj.time(obj.auTrigs(2)) obj.time(obj.auTrigs(2))], [-2 2], 'k--', 'LineWidth', lineThick)
        box off
        axis([auTimeRange min(obj.rawMic) max(obj.rawMic)])
        title('Raw Audapter Microphone')

        xlabel('Audapter Time (s)')
        set(gca,'FontName', 'Arial',...
                'FontSize', 14,...
                'FontWeight','bold')

        % Raw Audapter Headphones
        axes(ha(3))
        plot(obj.time, obj.rawHead)
        hold on
        plot([obj.voiceOnsetT obj.voiceOnsetT], [-2 2], 'Color', voiceOnsetColor, 'LineStyle', '--', 'LineWidth', lineThick)
        hold on
        plot([obj.time(obj.auTrigs(1)) obj.time(obj.auTrigs(1))], [-2 2], 'k--', 'LineWidth', lineThick)
        hold on
        plot([obj.time(obj.auTrigs(2)) obj.time(obj.auTrigs(2))], [-2 2], 'k--', 'LineWidth', lineThick)
        box off
        axis([auTimeRange min(obj.rawMic) max(obj.rawMic)])
        title('Raw Audapter Headphones')

        xlabel('Audapter Time (s)')
        set(gca,'FontName', 'Arial',...
                'FontSize', 14,...
                'FontWeight','bold')

        annotation('textbox', [0.8 0.9 0.1 0.1], 'String', {['MiHe Delay: ' num2str(obj.AuMHDelay) 's'];...
                                                            ['AuNI Delay: ' num2str(obj.AuNIDelay) 's']},...
                                                 'EdgeColor', 'none',...
                                                 'FontSize', 14,...
                                                 'FontWeight', 'Bold');

        end
    end
end

