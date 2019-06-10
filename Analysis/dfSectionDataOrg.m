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
        
        timeSec
        sigsSec
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
            
            % Time vector correspnding to the sectioned signals
            obj.timeSec = linspace(obj.preEveT, obj.posEveT, obj.numSampSec);
        end
        
        function obj = sectionData(obj)

            OnsetSecs  = [];
            OffsetSecs = [];
            if obj.numTrial > 0
                for ii = 1:obj.numTrial
                    OnsetT   = obj.trigs(ii, 1); % Onset time
                    OffsetT  = obj.trigs(ii, 2); % Offset time

                    OnsetTSt = round(OnsetT - obj.preEveT, 3);   % PreOnset time, rounded to nearest ms
                    OnsetTSp = round(OnsetT + obj.posEveT, 3);   % PostOnset time, rounded to nearest ms
                    OnsetSpan = obj.time >= OnsetTSt & obj.time <= OnsetTSp; % Indices corresponding to Onset period

                    OffsetTSt = round(OffsetT - preEve, 3); % PreOffset time, rounded to nearest ms
                    OffsetTSp = round(OffsetT + posEve, 3); % PostOffset time, rounded to nearest ms
                    OffsetSpan = obj.time >= OffsetTSt & obj.time <= OffsetTSp; % Indices corresponding to Offset period

                    OnsetSec  = obj.sigs(OnsetSpan, ii);  % Data sectioned around Onset
                    OffsetSec = obj.sigs(OffsetSpan, ii); % Data sectioned around Offset

                    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
                    OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
                end
            end

            obj.sigsSec(:,:,1) = OnsetSecs;  % 1st 3D layer
            obj.sigsSec(:,:,2) = OffsetSecs; % 2nd 3D layer
        end
    end
end

