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
            
            obj = obj.sectionData;
            obj = obj.meanData;
        end
        
        function obj = sectionData(obj)

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

                    OnsetSec  = obj.sigs(OnsetSpan, ii);  % Data sectioned around Onset
                    OffsetSec = obj.sigs(OffsetSpan, ii); % Data sectioned around Offset

                    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
                    OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
                end
            end

            obj.sigsSec(:,:,1) = OnsetSecs;  % 1st 3D layer
            obj.sigsSec(:,:,2) = OffsetSecs; % 2nd 3D layer
        end
        
        function obj = meanData(obj)
            % Some simple statistics on the sectioned audio data. 
            % meanAudio is a vector containing the following information
            % meanAudio(1) = mean Onset pitch contour
            % meanAudio(2) = 95% CI of the mean Onset Pitch Contour
            % meanAudio(3) = mean Offset pitch contour
            % meanAudio(4) = 95% CI of the mean Offset Pitch Contour

            OnsetSecs  = obj.sigsSec(:,:,1);
            OffsetSecs = obj.sigsSec(:,:,2);

            meanOnset  = nanmean(OnsetSecs, 2);  % across columns
            meanOffset = nanmean(OffsetSecs, 2); % across columns

            stdOnset   = nanstd(OnsetSecs, 0, 2);  % across columns
            stdOffset  = nanstd(OffsetSecs, 0, 2); % across columns

            SEMOnset   = stdOnset/sqrt(obj.numTrial);  % Standard Error
            SEMOffset  = stdOffset/sqrt(obj.numTrial); % Standard Error

            % NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
            % NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

            obj.sigsSecM = [meanOnset SEMOnset meanOffset SEMOffset];
        end
        
        function obj = identifyBaselineValues(obj)
            prePertT      = obj.timeSec <= 0;                      % timeSec is aligned for timeSec = 0 to be Onset of pert
            obj.sigsBase  = nanmean(obj.sigsSec(prePertT,:,1), 1); % Per-trial baseline value  
            obj.sigsBaseM = nanmean(obj.sigsBase);                 % Mean trial baseline value
        end
    end
end