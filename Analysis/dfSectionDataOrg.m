classdef dfSectionDataOrg
    % For the time being, lets assume two trigger points, onset and offset
    
    properties
        
        coder         % The person (if any) who took the measurements
        curSess       % Combined Participant and Run information
        dataType      % Type of data measured 
        dataUnit      % Units of the data measured
        iterationType % Level of analysis (cross-trial//cross-participant)
        
        time          % Time vector for the length of the full signals
        sigs          % Data vector for full length signals
        trigs         % Trigger points (time) when onset and offset occur
        fs            % sampling rate of the time/data vectors
        
        % Basic Sectioning Trial variables
        numTrial      % Number of trials involved
        preEveT       % Amount of time pre-trigger
        posEveT       % Amount of time post-trigger
        eveTLen       % Time length to section over from (-pre:pos)
        numSampSec    % Number of samples in a section
        pVec          % Vector of samples starting from first point
        
        sigsLims      % limits of the full signal
        
        % Sectioned trial versions
        timeSec       % time vector of the sectioned signal
        sigsSec       % signal sectioned around the trigger
        
        sigsBase      % baseline value of each trial (from raw signal)
        sigsBaseM     % mean baseline value
        
        sigsNorm      % full length signal normalized by baseline
        sigsNormSec   % sectioned signal normalized by baseline
        sigsNormSecSv % sectioned signal normalized with tossed trials removed
        
        sigsSecM      % mean sectioned signals
        sigsSecMLims  % limits of the mean sectioned signal
        
        sigsMeanFig      % figure for displaying the 
        sigsMeanFigTitle
        
        svIdx
        removedTrialTracker
        
        % Plotting Variables
        dataColor1
        dataColor2
        dataColor3
        
        f0TraceColor
        
        OnsetOffsetAxes
        legendCurves
        legendLabels
        LgdObj
        
        figTextName
        
        sigsDynamics
    end
    
    methods
        function obj = dfSectionDataOrg(time, sigs, trigs, fs, varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            if isempty(varargin)
                dataInfo.coder   = '';
                dataInfo.curSess = '';
                dataInfo.sigType = '';
                dataInfo.units   = '';
                dataInfo.itrType = '';
            else
                dataInfo = varargin{1};
            end
            
            obj.coder         = dataInfo.coder;
            obj.curSess       = dataInfo.curSess;
            obj.dataType      = dataInfo.sigType;
            obj.dataUnit      = dataInfo.units;
            obj.iterationType = dataInfo.itrType;
            
            obj.time  = time;
            obj.sigs  = sigs;
            obj.trigs = trigs;
            obj.fs    = fs;
            
            [~, obj.numTrial] = size(obj.sigs);
            
            obj.preEveT = -0.5; % 500ms before trigger
            obj.posEveT = 1.0;  % 1000ms after trigger
            obj.eveTLen = obj.posEveT - obj.preEveT;
            obj.numSampSec = obj.eveTLen*obj.fs+1;
            
            obj.pVec = linspace(0, obj.numSampSec-1, obj.numSampSec);
            
            % Time vector corresponding to the sectioned signals
            obj.timeSec = linspace(obj.preEveT, obj.posEveT, obj.numSampSec);
            
            obj.dataColor1 = [55,126,184]/255; % Cerulean Blu
            obj.dataColor2 = [77,175,74]/255;  % Leaf Green
            obj.dataColor3 = [231,41,138]/255; % Bright Magenta
            
            obj.f0TraceColor = [0 0 1];
            obj.figTextName = 'Arial';
            
            obj.legendCurves = [];
            obj.legendLabels = {};
            
            obj.svIdx = [];
            obj.removedTrialTracker = [];
            
            % Mean the trials
%             obj.sigsSecM = obj.meanData(obj.sigsNormSec);
            
            % Identify limits of the mean trials
%             obj.sigsSecMLims = obj.identifyBounds;
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
        
        function obj = identifyBaselineValues(obj, sigsSec)
            prePertT      = obj.timeSec <= 0;                  % timeSec is aligned for timeSec = 0 to be Onset of pert
            obj.sigsBase  = nanmean(sigsSec(prePertT,:,1), 1); % Per-trial baseline value  
            obj.sigsBaseM = nanmean(obj.sigsBase);             % Mean trial baseline value
        end
        
        function obj = convertCentsData(obj)
            % This will normalize the full length signal into another full
            % length signal by the equation for converting a f0 recording
            % from Hz to cents
            
            obj.sigsNorm = [];
            for ii = 1:obj.numTrial
                norm = 1200*log2(obj.sigs(:,ii)./obj.sigsBase(ii));
                obj.sigsNorm  = cat(2, obj.sigsNorm, norm);
            end
        end
        
        function sigsSmooth = applySmoothing(obj, sigs, len)
            
            sigsSmooth = [];
            for ii = 1:obj.numTrial
                smoothed   = smooth(sigs(:,ii), len);   % 10 sample length smoothing
                sigsSmooth = cat(2, sigsSmooth, smoothed);
            end
        end
        
        function sigsZM = applyZeroMean(obj, sigs, base)
            sigsZM = sigs - base;
        end
        
        function obj = qualityCheckData(obj, sigSec, varargin)
            % Identify trials (mic) where participant had vocal fry
            % aka: f0 pitch miscalculation
            
            if isempty(varargin)
                curSavedIdxList = 1:obj.numtrial;
            else
                curSavedIdxList = varargin{1};
            end
            
            for ii = 1:obj.numTrial
                curSavedIdx = curSavedIdxList(ii);
                sigSecTrial = sigSec(:, ii, :);

                % Any points where pitch is > 500 // < -500 cents?
                MisCalc = find(sigSecTrial >= 500 | sigSecTrial <=  -500);

                if ~isempty(MisCalc)
%                     if obj.trialTypeSvt(ii) == 0
%                         type = 'Control';
%                     else
%                         type = 'Perturbed';
%                     end       
%                     fprintf('%s Trial %d (%s) excluded due to Miscalculated Pitch Trace\n', obj.curSess, svIdc, type)

                    removedTrial = {['Trial ' num2str(curSavedIdx)], 'Miscalculated pitch Trace'};
                    obj.removedTrialTracker = cat(1, obj.removedTrialTracker, removedTrial);
                else
                    obj.svIdx = cat(1, obj.svIdx, ii); % Keep the index of saved trials (not removed)
                end
            end
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

            sigsSecM.ON.mean = nanmean(OnsetSecs, 2);  % across columns
            sigsSecM.OF.mean = nanmean(OffsetSecs, 2); % across columns

            sigsSecM.ON.STD = nanstd(OnsetSecs, 0, 2);  % across columns
            sigsSecM.OF.STD = nanstd(OffsetSecs, 0, 2); % across columns

            sigsSecM.ON.SEM = sigsSecM.ON.STD/sqrt(obj.numTrial);  % Standard Error
            sigsSecM.OF.SEM = sigsSecM.OF.STD/sqrt(obj.numTrial); % Standard Error

            sigsSecM.ON.NCI = 1.96*sigsSecM.ON.SEM;  % 95% Confidence Interval
            sigsSecM.OF.NCI = 1.96*sigsSecM.OF.SEM; % 95% Confidence Interval
        end
        
        function obj = identifyBounds(obj)
            
            % Full Trial Limits
            uLSigs = max(obj.sigs);
            lLSigs = min(obj.sigs);
            
            uL = round(max(uLSigs)) + 10;
            lL = round(min(lLSigs)) - 10;
            
            obj.sigsLims = [obj.time(1) obj.time(end) lL uL];
            
            % Section Mean Limits
            [lwBoundOn, upBoundOn] = obj.MaxMinBounds(obj.sigsSecM.ON.mean, obj.sigsSecM.ON.NCI, 10);
            [lwBoundOf, upBoundOf] = obj.MaxMinBounds(obj.sigsSecM.OF.mean, obj.sigsSecM.OF.NCI, 10);

            lwBoundM = obj.CompareAndChooseBounds(lwBoundOn, lwBoundOf, 'min');
            upBoundM = obj.CompareAndChooseBounds(upBoundOn, upBoundOf, 'max');

            obj.sigsSecMLims = [obj.timeSec(1) obj.timeSec(end) lwBoundM upBoundM];
        end
        
        function [lwBound, upBound] = MaxMinBounds(obj, audio, audioErr, buff)
            % MaxMinBounds takes an time-series measure and the measured
            % error and identifies the max and min values for the 
            % time-series trace. It identifies how wide the bounds needs to
            % be fully show the trace with some buffer

            [~, Imax] = max(audio);
            upBound   = round(audio(Imax) + audioErr(Imax) + buff);
            [~, Imin] = min(audio);
            lwBound   = round(audio(Imin) - audioErr(Imin) - buff);
        end

        function bound = CompareAndChooseBounds(obj, bound1, bound2, type)
            % CompareAndChooseBounds does a simple and quick check between two proposed
            % bounds and decides which is greater/lesser. Since I am often comparing
            % two sets of things, this helps to keep data ranges consistent and
            % symmetrical. 

            allBounds = [bound1 bound2];

            if strcmp(type, 'max')
                bound = max(allBounds);
            elseif strcmp(type, 'min')
                bound = min(allBounds);
            else
                bound = allBounds(1);
            end

        end
        
        function obj = drawSigsSecM(obj, varargin)
            % drawSigsSecM(obj) plots onset and offset sectioned signals
            % manipulated in this class.

            if isempty(varargin)
                plotStimWindows = 0;
            else
                plotStimWindows = 1;
                stimWindowProp = varargin{1};
            end
            
            % Figure properties
            plotpos = [10 100];
            plotdim = [1600 600];
            OnsetOffsetMeanDataFig = figure('Color', [1 1 1]);
            set(OnsetOffsetMeanDataFig, 'Position', [plotpos plotdim], 'PaperPositionMode','auto')

%             AD = res.audioDynamics;
%             if ~isempty(AD)
%                 statSM = round(AD.respVarM(2), 1);
%                 statRM = round(AD.respVarM(3), 1);
%                 statRP = round(AD.respVarM(4));
%             else
%                 statSM = '';
%                 statRM = '';
%                 statRP = '';
%             end
%             curSess(strfind(curSess, '_')) = ' ';

            % Trigger Line properties
            trigLineX   = [0 0];
            trigLineY   = [-10000 10000];
            trigLineC   = [0.3 0.3 0.3];
            
            % Zero Line properties
            zeroLineX = [obj.timeSec(1) obj.timeSec(end)];
            zeroLineY = [0 0];
            zeroLineC = [0.3 0.3 0.3];
            
            % Axis properties
            titleFontSize  = 14;
            axisLSize      = 14;
            lineThick      = 4;

            % Subplot dimensions and spacing
            obj.OnsetOffsetAxes = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

            % Onset of Perturbation
            axes(obj.OnsetOffsetAxes(1))
            if plotStimWindows == 1
                onsetStimWinC = [0 1 0];
                onsetStimWinAx = [stimWindowProp.meanOnsetLag stimWindowProp.meanOnsetLag+stimWindowProp.meanOnsetRise];
                onsetStimWinAy = [600 600];
                area(onsetStimWinAx, onsetStimWinAy, -600, 'FaceColor', onsetStimWinC, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);
                hold on
            end
                        
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            dHOn = shadedErrorBar(obj.timeSec, obj.sigsSecM.ON.mean, obj.sigsSecM.ON.NCI, 'lineprops', {'color', obj.dataColor1}, 'transparent', 1);

            set(dHOn.mainLine, 'LineWidth', lineThick)
            xlabel('Time (s)', 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold')
            title('Onset of Perturbation', 'FontName', obj.figTextName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off

            set(gca,'FontName', obj.figTextName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold',...
                    'LineWidth', 2)
            hold off
            
            % Offset of Perturbation
            axes(obj.OnsetOffsetAxes(2))
            
            if plotStimWindows == 1
                offsetStimWinC = [1 0 0];
                offsetStimWinAx = [stimWindowProp.meanOffsetLag stimWindowProp.meanOffsetLag+stimWindowProp.meanOffsetRise];
                offsetStimWinAy = [600 600];
                area(offsetStimWinAx, offsetStimWinAy, -600, 'FaceColor', offsetStimWinC, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);
                hold on
            end
            
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            dHOf = shadedErrorBar(obj.timeSec, obj.sigsSecM.OF.mean, obj.sigsSecM.OF.NCI, 'lineprops', {'color', obj.dataColor1}, 'transparent', 1);
              
            obj.legendCurves = cat(2, obj.legendCurves, dHOf.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, ['Line 1: ' num2str(obj.numTrial) ' ' obj.iterationType]);
            
            set(dHOf.mainLine, 'LineWidth', lineThick)
            xlabel('Time (s)', 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold')
            title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off
            hold off
            set(gca,'FontName', obj.figTextName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold',...
                    'LineWidth', 2,...
                    'YAxisLocation', 'right');

            sup = suptitle({obj.curSess});
            set(sup, 'FontName', obj.figTextName,...
                     'FontSize', titleFontSize,...
                     'FontWeight','bold')
                 
            annotation('textbox',[0.88 0.88 0.45 0.1],...
                       'string', ['Coder: ' obj.coder],...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',18,...
                       'FontName','Arial');

            obj.LgdObj = legend(obj.legendCurves, obj.legendLabels,...
                                'Box', 'off',...
                                'Edgecolor', [1 1 1],...
                                'FontSize', 12,...
                                'FontWeight', 'bold',...
                                'Position', [0.7 0.91 0.05 0.05]);         


            obj.sigsMeanFig      = OnsetOffsetMeanDataFig;
            obj.sigsMeanFigTitle = [obj.curSess '_InterTrialMean' obj.coder '.jpg'];
        end
        
        function obj = drawSigsSecM_Onset(obj, shadeFlag, varargin)
            % drawSigsSecM_Onset(obj) plots onset sectioned signals

            if isempty(varargin)
                plotStimWindows = 0;
            else
                plotStimWindows = 1;
                stimTraceProp = varargin{1};
            end
            
            % Figure properties
            plotpos = [10 100];
            plotdim = [800 600];
            OnsetMeanDataFig = figure('Color', [1 1 1]);
            set(OnsetMeanDataFig, 'Position', [plotpos plotdim], 'PaperPositionMode','auto')
            
            set(OnsetMeanDataFig,'defaultAxesColorOrder',[obj.dataColor1; obj.dataColor3]);

            % Trigger Line properties
            trigLineX   = [0 0];
            trigLineY   = [-10000 10000];
            trigLineC   = [0.3 0.3 0.3];
            
            % Zero Line properties
            zeroLineX = [obj.timeSec(1) obj.timeSec(end)];
            zeroLineY = [0 0];
            zeroLineC = [0 0 0];
            
            % Pressure Line properties
            pressureC = [255 87 51]/255; %Auburn
            
            % Axis properties
            axisLSize      = 20;
            lineThick      = 4;

            % Onset of Perturbation
            if plotStimWindows == 1
                minPres = min(stimTraceProp.meanPres.OF.mean); maxPres = max(stimTraceProp.meanPres.ON.mean);
                m = (obj.sigsSecMLims(4) - obj.sigsSecMLims(3))/(maxPres - minPres);
                b = obj.sigsSecMLims(3) - m*minPres;
                presTrace = plot(stimTraceProp.timeP, m*stimTraceProp.meanPres.ON.mean+b, 'color', pressureC, 'LineStyle', ':', 'LineWidth', lineThick);
                obj.legendCurves = cat(2, obj.legendCurves, presTrace);
                obj.legendLabels = cat(2, obj.legendLabels, 'Balloon pressure');
                hold on
            end
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            
            % Draw the Onset Trace
            if shadeFlag == 1
                dHOn = shadedErrorBar(obj.timeSec, obj.sigsSecM.ON.mean, obj.sigsSecM.ON.NCI, 'lineprops', {'color', obj.dataColor1, 'LineWidth', lineThick}, 'transparent', 1);
                obj.legendCurves = cat(2, obj.legendCurves, dHOn.mainLine);
                obj.legendLabels = cat(2, obj.legendLabels, 'Laryngeal movement index');
            else
                dHOn = plot(obj.timeSec, obj.sigsSecM.ON.mean, 'color', 'k', 'LineWidth', lineThick);
                obj.legendCurves = cat(2, obj.legendCurves, dHOn);
                obj.legendLabels = cat(2, obj.legendLabels, 'Laryngeal movement index');
            end
            
            xlabel('Time (s)', 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold', 'Color', obj.dataColor1)
%             title('Onset of Perturbation', 'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off

            set(gca,'FontName', obj.figTextName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold',...
                    'LineWidth', 2,...
                    'YColor', obj.dataColor1)
            hold off
            
            obj.LgdObj = legend(obj.legendCurves, obj.legendLabels,...
                    'FontName', obj.figTextName,...
                    'FontSize', 16,...
                    'FontWeight', 'bold',...
                    'Position', [0.608 0.145 0.1 0.1],...
                    'Box', 'Off'); 

            obj.sigsMeanFig      = OnsetMeanDataFig;
            obj.sigsMeanFigTitle = [obj.curSess '_InterTrialMean' obj.coder '.jpg'];
        end
        
        function obj = drawSigsSecM_Offset(obj, shadeFlag, varargin)
            % drawSigsSecM_Onset(obj) plots onset sectioned signals

            if isempty(varargin)
                plotStimWindows = 0;
            else
                plotStimWindows = 1;
                stimTraceProp = varargin{1};
            end
            
            % Figure properties
            plotpos = [10 100];
            plotdim = [800 600];
            OffsetMeanDataFig = figure('Color', [1 1 1]);
            set(OffsetMeanDataFig, 'Position', [plotpos plotdim], 'PaperPositionMode','auto')
            
            set(OffsetMeanDataFig,'defaultAxesColorOrder',[obj.dataColor1; obj.dataColor3]);

            % Trigger Line properties
            trigLineX   = [0 0];
            trigLineY   = [-10000 10000];
            trigLineC   = [0.3 0.3 0.3];
            
            % Zero Line properties
            zeroLineX = [obj.timeSec(1) obj.timeSec(end)];
            zeroLineY = [0 0];
            zeroLineC = [0 0 0];
        
            % Pressure Line properties
            pressureC = [255 87 51]/255; %Auburn
            
            % Axis properties
            axisLSize      = 20;
            lineThick      = 4;

            % Onset of Perturbation
            if plotStimWindows == 1
                minPres = min(stimTraceProp.meanPres.OF.mean); maxPres = max(stimTraceProp.meanPres.ON.mean);
                m = (obj.sigsSecMLims(4) - obj.sigsSecMLims(3))/(maxPres - minPres);
                b = obj.sigsSecMLims(3) - m*minPres;
                presTrace = plot(stimTraceProp.timeP, m*stimTraceProp.meanPres.OF.mean+b, 'color', pressureC, 'LineStyle', ':', 'LineWidth', lineThick);
                obj.legendCurves = cat(2, obj.legendCurves, presTrace);
                obj.legendLabels = cat(2, obj.legendLabels, 'Balloon pressure');
                hold on
            end
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            
            % Draw the Onset Trace
            if shadeFlag == 1
                dHOn = shadedErrorBar(obj.timeSec, obj.sigsSecM.OF.mean, obj.sigsSecM.OF.NCI, 'lineprops', {'color', obj.dataColor1, 'LineWidth', lineThick}, 'transparent', 1);
                obj.legendCurves = cat(2, obj.legendCurves, dHOn.mainLine);
                obj.legendLabels = cat(2, obj.legendLabels, 'Laryngeal movement index');
            else
                dHOn = plot(obj.timeSec, obj.sigsSecM.OF.mean, 'color', 'k', 'LineWidth', lineThick);
                obj.legendCurves = cat(2, obj.legendCurves, dHOn);
                obj.legendLabels = cat(2, obj.legendLabels, 'Laryngeal movement index');
            end
            
            xlabel('Time (s)', 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', obj.figTextName, 'FontSize', axisLSize, 'FontWeight', 'bold', 'Color', obj.dataColor1)
%             title('Onset of Perturbation', 'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off

            set(gca,'FontName', obj.figTextName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold',...
                    'LineWidth', 2,...
                    'YColor', obj.dataColor1)
            hold off
            
            obj.LgdObj = legend(obj.legendCurves, obj.legendLabels,...
                    'FontName', obj.figTextName,...
                    'FontSize', 12,...
                    'FontWeight', 'bold',...
                    'Position', [0.635 0.165 0.1 0.1]); 

            obj.sigsMeanFig      = OffsetMeanDataFig;
            obj.sigsMeanFigTitle = [obj.curSess '_InterTrialMean' obj.coder '.jpg'];
        end
        
        function obj = appendFigure(obj, sigsSecM, flag)
            % obj = appendFigure(obj, sigsSecM, flag) allows you to add
            % more meaned section traces to the figure which has already
            % been created for the defauly trace.
            
            if flag == 2
                color = obj.dataColor2;
            else
                color = obj.dataColor3;
            end
            
            axes(obj.OnsetOffsetAxes(1))
            hold on
            onH = shadedErrorBar(obj.timeSec, sigsSecM.ON.mean, sigsSecM.ON.NCI, 'lineprops', {'color', color}, 'transparent', 1);
            set(onH.mainLine, 'LineWidth', 4)
            hold off
            
            axes(obj.OnsetOffsetAxes(2))
            hold on
            ofH = shadedErrorBar(obj.timeSec, sigsSecM.OF.mean, sigsSecM.OF.NCI, 'lineprops', {'color', color}, 'transparent', 1);
            set(ofH.mainLine, 'LineWidth', 4)
            hold off
            
            obj.legendCurves = cat(2, obj.legendCurves, ofH.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, ['Line ' num2str(flag) ': ' num2str(obj.numTrial) ' ' obj.iterationType]);
            
            legend(obj.legendCurves, obj.legendLabels);
        end
        
        function saveSigsSecMFig(obj, plotFolder)
            
            saveFileName = fullfile(plotFolder, obj.sigsMeanFigTitle);
            export_fig(saveFileName, '-r300')
        end
        
        function ir = initInflationResponseStruct(obj)

            ir.time     = [];
            ir.onset    = [];

            ir.iAtOnset = []; % Index where t = 0
            ir.tAtOnset = []; % Time at t = 0
            ir.vAtOnset = []; % f0 value at t = 0

            ir.iPostOnsetR = []; % Range of indices between t = 0ms and t = 200ms;
            ir.iAtMin      = []; % Index at min f0 value in PostOnsetR
            ir.tAtMin      = []; % Time at min f0 value in PostOnsetR
            ir.vAtMin      = []; % Min f0 value in PostOnsetR
            ir.stimMag     = []; % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

            ir.iAtResp = []; % Index of f0 value when participant 'fully' responded...right now = last value in section
            ir.tAtResp = []; % Time at f0 value when participant 'fully' responded
            ir.vAtResp = []; % f0 value when participant 'fully' responded 
            ir.respMag = []; % vAtResp - vAtMin   ...distance traveled
            ir.respPer = []; % Percent change from stimMag to respMag
        end
        
        function audioDynamics_Somato = InflationResponse(obj)
        % [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
        % Identifies the relevant pitch contour characteristics that are important
        % for deciding how a participant responded to the inflation of the balloon
        % during production. iR is a structure representing the result variables
        % from studying the inflation response (iR). The prefix letter denotes
        % whether the variable is a index (i), a time (t), or a value (v). 
        %
        % secTime:  Vector of time points corresponding to the sectioned data (numSamp)
        % secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
        %           The 1st 3D layer are Onset Sections
        %           The 2nd 3D later are Offset Sections
        %
        % respVar: Matrix of per trial iR results (numTrial x 4)
        %          respVar(:,1) = Time of the minimum f0 value in the sec
        %          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
        %          respVar(:,3) = Value of f0 at end of sec (response magnitude)
        %          respVar(:,4) = ABS percent change of stim and response
        %          magnitude (response percentage)
        % respVarM:    Vector of mean trial values from respVarm (1x4)
        % respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
        % InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

        ir = obj.initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
        ir.time     = obj.timeSec;          % Time Interval for the sectioned trials (-0.5->1.0s)
        ir.onset    = obj.sigsSecM.ON.mean;   % f0 Trace sectioned around pert Onset.

        ir.iAtOnset = find(ir.time == 0);
        ir.tAtOnset = 0;                     % duh
        ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

        ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
        [minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

        % StimMag
        ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
        ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
        ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
        ir.stimMag = abs(ir.vAtMin - ir.vAtOnset); % Distance traveled from onset to min value

        % RespMag
        ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
        ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
        ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
        ir.tAtResp      = mean(ir.tAtRespRange);
        ir.vAtResp      = mean(ir.vAtRespRange);
        ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

        % RespPer
        ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

        % Add to the audioDynamics struct
        respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
        audioDynamics_Somato.respVarM = respVarM;
        end
        
        function appendFigureDynamics_Onset(obj, time, sec)
            yyaxis right
            f0Trace = shadedErrorBar(time, sec.ON.mean, sec.ON.NCI, 'lineprops', {'color', obj.f0TraceColor, 'LineWidth', 4}, 'transparent', 1);
            obj.legendCurves = cat(2, obj.legendCurves, f0Trace.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, '{\it f}_o');
            obj.LgdObj = legend(obj.legendCurves, obj.legendLabels);
            ylabel('{\it f}_o (cents)')
            axis([-0.5 1.0 -108 158])
            set(gca,'YColor', obj.f0TraceColor)
            yyaxis left
            axis([-0.5 1.0 -10.8 15.8])
            set(gca,'YColor', obj.dataColor1)
        end
        
        function appendFigureDynamics_Offset(obj, time, sec)
            yyaxis right
            f0Trace = shadedErrorBar(time, sec.OF.mean, sec.OF.NCI, 'lineprops', {'color', obj.f0TraceColor, 'LineWidth', 4}, 'transparent', 1);
            obj.legendCurves = cat(2, obj.legendCurves, f0Trace.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, '{\it f}_o');
            obj.LgdObj = legend(obj.legendCurves, obj.legendLabels);
            ylabel('{\it f}_o(cents)')
            axis([-0.5 1.0 -108 158])
            set(gca,'YColor', obj.f0TraceColor)
            yyaxis left
            axis([-0.5 1.0 -10.8 15.8])
            set(gca,'YColor', obj.dataColor1)
        end
    end
end