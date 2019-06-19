classdef MicHeadAlignProcess
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
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
        rawMicDS    % Raw microphone signal, downsampled
        
        frameDel
        rmsThresh
        
        rawMicNI
        fsNI
        trialLenNI
        trialTimeNI
        tNI
        
        numSamp
        
        rmsVoiceFrame
        voiceOnsetInd
        voiceOnsetT
        
        preVOnsetRMS
        
        AuNIDelay
        AuNIDelayP
        adjustedDelay
        
        AuMHDelay
        AuMHDelayP
        
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
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
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
            obj.rawMicNI    = trialVar.rawMicNI;
            obj.fsNI        = analysisVar.sRateNi;
            obj.trialLenNI  = length(obj.rawMicNI);
            obj.trialTimeNI = obj.trialLenNI/obj.fsNI;
            obj.tNI         = linspace(0, obj.trialTimeNI, obj.trialLenNI);
            
            obj.numSamp   = obj.trialTimeNI*obj.fs;

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
            obj.AuNIDelayP = obj.AuNIDelay*obj.fs;                          % Convert to points

            % Adjust Triggers against NIDAQ if a Laryngeal Pert Exp.
            % Otherwise adjust based on VoiceOnset, which happens during the PSR
            if strcmp(obj.expType(1:3), 'Som')
                obj.adjustedDelay = obj.AuNIDelayP;
            else
                obj.adjustedDelay = obj.voiceOnsetInd;
            end
            obj.auTrigsAuNi = obj.auTrigs + obj.adjustedDelay;
            obj.auTimesAuNi = obj.time(obj.auTrigsAuNi);

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
        
        function drawPreProcessDiagnostic(obj, check)

        plotPos = [10 10];
        plotDim = [1200 900];

        MHFig = figure('Color', [1 1 1]);
        set(MHFig, 'Position', [plotPos plotDim],'PaperPositionMode','auto')

        ha = tight_subplot(2,1,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

        axes(ha(1))
        plot(obj.time, obj.rawMic)
        hold on
        plot(obj.time, obj.env, 'y')
        hold on
        plot([obj.voiceOnsetT obj.voiceOnsetT], [-0.2 0.2])
        hold on
        plot([obj.time(obj.auTrigs(1)) obj.time(obj.auTrigs(1))], [-0.2 0.2], 'k--')
        hold on
        plot([obj.time(obj.auTrigs(2)) obj.time(obj.auTrigs(2))], [-0.2 0.2], 'k--')
        box off
        axis([0 6 -0.25 0.25])

        if check == 1
            axes(ha(2))
            plot(obj.time, obj.rawHead)
            hold on
            plot([obj.voiceOnsetT obj.voiceOnsetT ], [-0.2 0.2])
            hold on
            plot([obj.time(obj.auTrigs(1)) obj.time(obj.auTrigs(1))], [-0.2 0.2], 'k--')
            hold on
            plot([obj.time(obj.auTrigs(2)) obj.time(obj.auTrigs(2))], [-0.2 0.2], 'k--')
            box off
            axis([0 6 -0.25 0.25])
            title(num2str(obj.AuMHDelay))
        else
            axes(ha(2))
            plot(obj.tNI, obj.rawMicNI)
            title(num2str(obj.AuNIDelay))
            axis([0 4 -0.25 0.25])
            box off
        end
        end
    end
end
