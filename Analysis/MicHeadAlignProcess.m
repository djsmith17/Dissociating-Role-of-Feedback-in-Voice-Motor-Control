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
        t           % time vector of the recording
        auTrigs     % trigger points from Audapter recording
        auTrigsAuNi % trigger points from Audapter recording aligned with NIDAQ triggers
        rawMicDS    % Raw microphone signal, downsampled
        
        frameDel
        rmsThresh
        
        rawMicNI
        fsNI
        trialLenNI
        trialTimeNI
        tNI
        
        numSamp
        envThresh
        
        env
        maxPeak
        threshIdx
        rmsVoiceInd
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
            
            obj.expType   = analysisVar.expType;
            obj.AudFB     = trialVar.AudFB;
            obj.rawMic    = double(trialVar.rawMic);
            obj.rawHead   = double(trialVar.rawHead);
            obj.rms       = double(trialVar.rms);
            obj.fs        = analysisVar.sRate;
            obj.frameLen  = analysisVar.frameLen;
            obj.trialLen  = length(obj.rawMic);
            obj.trialTime = obj.trialLen/obj.fs;
            obj.t         = linspace(0, obj.trialTime, obj.trialLen);
            obj.auTrigs   = trialVar.auTrigs;
            obj.frameDel  = 7;
            obj.rmsThresh = 0.011;
            
            % NIDAQ Recording Specifics
            obj.rawMicNI    = trialVar.rawMicNI;
            obj.fsNI        = analysisVar.sRateNi;
            obj.trialLenNI  = length(obj.rawMicNI);
            obj.trialTimeNI = obj.trialLenNI/obj.fsNI;
            obj.tNI         = linspace(0, obj.trialTimeNI, obj.trialLenNI);
            
            obj.numSamp   = obj.trialTimeNI*obj.fs;
            obj.envThresh = 0.30;
            
            % Find the envelope of the microphone signal
            obj.env = obj.calcEnvelope(obj.rawMic, obj.fs);
            
            % Largest peak in the envelope theoretically occurs during voicing
            obj.maxPeak = max(obj.env);

            % Find values that are within threshold of max 'voicing' value
            obj.threshIdx = find(obj.env > obj.envThresh*obj.maxPeak); 
            
            % First index of the theoretical useable signal (Voice onset)
            obj.rmsVoiceInd   = find(obj.rms > obj.rmsThresh);
            obj.voiceOnsetInd = (obj.rmsVoiceInd(1) - obj.frameDel)*obj.frameLen;
            if obj.voiceOnsetInd <= 0 % If they started speaking IMMEDIATELY...can't have index of 0
                obj.voiceOnsetInd = 1;
            end
            obj.voiceOnsetT   = obj.t(obj.voiceOnsetInd);
            
            % Evaluate the pre-voice onset rms (used to identify voice breaks)
            obj.preVOnsetRMS = evalPreVoiceRMS(obj);
            
            % Find the delay between NIDAQ recording and Audapter recording
            obj.rawMicDS   = resample(obj.rawMic, obj.fsNI, obj.fs);         % Downsample the Audapter recording
            obj.AuNIDelay  = obj.xCorrTimeLag(obj.rawMicDS, obj.rawMicNI, obj.fsNI); % Perform xCorr between NIDAQ and Audapter. Expect that NIDAQ leads Audapter
            obj.AuNIDelayP = obj.AuNIDelay*obj.fs;                          % Convert to points

            % Adjust Triggers against NIDAQ only if we are using Laryngeal Pert Exp.
            % Otherwise adjsut based on VoiceOnset, which is what Audapter does in PSR
            if strcmp(obj.expType(1:3), 'Som')
                obj.adjustedDelay = obj.AuNIDelayP;
                obj.auTrigsAuNi = obj.auTrigs + obj.adjustedDelay;
            else
                obj.adjustedDelay = obj.voiceOnsetInd;
                obj.auTrigsAuNi = obj.auTrigs + obj.adjustedDelay;
            end

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
        
        function env = calcEnvelope(obj, audio, fs)
            cutoffF = 40;
            [B, A] = butter(4, cutoffF/(fs/2));

            % Envelope the signal by low-pass filtering (change in amplitude/time ~RMS)
            env = filter(B, A, abs(audio)); 
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

            analysisPerFO  = obj.rms(obj.analysisFrames) < obj.preVOnsetRMS;
            postVOnsetFO   = obj.rms(postVOnsetFrames) < obj.preVOnsetRMS;

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
        plot(obj.t, obj.rawMic)
        hold on
        plot(obj.t, obj.env, 'y')
        hold on
        plot([obj.voiceOnsetT obj.voiceOnsetT ], [-0.2 0.2])
        hold on
        plot([obj.t(obj.auTrigs(1)) obj.t(obj.auTrigs(1))], [-0.2 0.2], 'k--')
        hold on
        plot([obj.t(obj.auTrigs(2)) obj.t(obj.auTrigs(2))], [-0.2 0.2], 'k--')
        box off
        axis([0 6 -0.25 0.25])

        if check == 1
            axes(ha(2))
            plot(obj.t, obj.rawHead)
            hold on
            plot([obj.voiceOnsetT obj.voiceOnsetT ], [-0.2 0.2])
            hold on
            plot([obj.t(obj.auTrigs(1)) obj.t(obj.auTrigs(1))], [-0.2 0.2], 'k--')
            hold on
            plot([obj.t(obj.auTrigs(2)) obj.t(obj.auTrigs(2))], [-0.2 0.2], 'k--')
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

