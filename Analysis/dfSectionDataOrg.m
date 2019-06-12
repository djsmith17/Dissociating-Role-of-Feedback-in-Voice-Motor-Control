classdef dfSectionDataOrg
    % For the time being, lets assume two trigger points, onset and offset
    
    properties
        
        time
        sigs
        trigs
        fs
        
        numTrial
        preEveT
        posEveT
        eveTLen
        numSampSec
        pVec
        
        timeSec
        sigsSec
        sigsSecM
        
        sigsBase
        sigsBaseM
        
        sigsNorm
        sigsNormSec
    end
    
    methods
        function obj = dfSectionDataOrg(time, sigs, trigs, fs)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.time  = time;
            obj.sigs  = sigs;
            obj.trigs = trigs;
            obj.fs    = fs;
            
            [~, obj.numTrial] = size(obj.sigs);
            
            obj.preEveT = -0.5; % 500ms before trigger
            obj.posEveT = 1.0;  % 1000ms after trigger
            obj.eveTLen = obj.posEveT - obj.preEveT;
            obj.numSampSec = obj.eveTLen*obj.fs;
            
            obj.pVec = linspace(0, obj.numSampSec-1, obj.numSampSec);
            
            % Time vector correspnding to the sectioned signals
            obj.timeSec = linspace(obj.preEveT, obj.posEveT, obj.numSampSec);
            
            % Section raw f0 around onset and offset
            obj.sigsSec  = obj.sectionData(obj.sigs);
            
            % Identify baseline values 
            
            obj.sigsSecM = obj.meanData(obj.sigsSec);
        end
        
        function sigsSec = sectionData(obj, sigs)

            OnsetSecs  = [];
            OffsetSecs = [];
            if obj.numTrial > 0
                for ii = 1:obj.numTrial
                    OnsetT   = obj.trigs(ii, 1); % Onset time
                    OffsetT  = obj.trigs(ii, 2); % Offset time

                    OnsetTSt = round(OnsetT + obj.preEveT, 3);   % PreOnset time, rounded to nearest ms
                    OnsetTStLeast = find(obj.time <= OnsetTSt);                  
                    OnsetSpan = OnsetTStLeast(end) + obj.pVec; % Indices corresponding to Onset period

                    OffsetTSt = round(OffsetT + obj.preEveT, 3); % PreOffset time, rounded to nearest ms
                    OffsetTStLeast = find(obj.time <= OffsetTSt);                  
                    OffsetSpan = OffsetTStLeast(end) + obj.pVec; % Indices corresponding to Onset peri

                    OnsetSec  = sigs(OnsetSpan, ii);  % Data sectioned around Onset
                    OffsetSec = sigs(OffsetSpan, ii); % Data sectioned around Offset

                    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
                    OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
                end
            end

            sigsSec(:,:,1) = OnsetSecs;  % 1st 3D layer
            sigsSec(:,:,2) = OffsetSecs; % 2nd 3D layer
        end
        
        function sigsSecM = meanData(obj, sigsSec)
            % Some simple statistics on the sectioned audio data. 
            % meanAudio is a vector containing the following information
            % meanAudio(1) = mean Onset pitch contour
            % meanAudio(2) = 95% CI of the mean Onset Pitch Contour
            % meanAudio(3) = mean Offset pitch contour
            % meanAudio(4) = 95% CI of the mean Offset Pitch Contour

            OnsetSecs  = sigsSec(:,:,1);
            OffsetSecs = sigsSec(:,:,2);

            meanOnset  = nanmean(OnsetSecs, 2);  % across columns
            meanOffset = nanmean(OffsetSecs, 2); % across columns

            stdOnset   = nanstd(OnsetSecs, 0, 2);  % across columns
            stdOffset  = nanstd(OffsetSecs, 0, 2); % across columns

            SEMOnset   = stdOnset/sqrt(obj.numTrial);  % Standard Error
            SEMOffset  = stdOffset/sqrt(obj.numTrial); % Standard Error

            % NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
            % NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

            sigsSecM = [meanOnset SEMOnset meanOffset SEMOffset];
        end
        
        function obj = identifyBaselineValues(obj, sigsSec)
            prePertT      = obj.timeSec <= 0;                      % timeSec is aligned for timeSec = 0 to be Onset of pert
            obj.sigsBase  = nanmean(sigsSec(prePertT,:,1), 1); % Per-trial baseline value  
            obj.sigsBaseM = nanmean(obj.sigsBase);                 % Mean trial baseline value
        end
        
        function obj = convertCentsData(obj)
            obj.sigsNorm = [];
            for ii = 1:obj.numTrial
                norm = 1200*log2(obj.sigs(:,ii)./obj.sigsBase(ii));
                obj.sigsNorm  = cat(2, obj.sigsNorm, norm);
            end
        end
        
        function obj = qualityCheckData(obj)
            % Identify trials (mic) where participant had vocal fry
            % aka: f0 pitch miscalculation
            
            for ii = 1:obj.numTrial

                mic     = obj.audioMf0_norm(timeInd, ii);

                % objy points where pitch is > 500 // < -500 cents?
                MisCalcf0 = find(mic >= 500 | mic <=  -500);

                if ~isempty(MisCalcf0)
                    if obj.trialTypeSvt(ii) == 0
                        type = 'Control';
                    else
                        type = 'Perturbed';
                    end       
                    fprintf('%s Trial %d (%s) excluded due to Miscalculated Pitch Trace\n', obj.curSess, svIdc, type)

                    removedTrial = {['Trial ' num2str(svIdc)], 'Miscalculated pitch Trace'};
                    obj.removedTrialTracker = cat(1, obj.removedTrialTracker, removedTrial);
                else
                    obj.subSvIdx = cat(1, obj.subSvIdx, ii); % Keep the index of the trials that are not removed
                end
            end
        end
    end
end