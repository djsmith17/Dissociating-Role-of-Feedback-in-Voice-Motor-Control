classdef dfSectionDataOrg
    % For the time being, lets assume two trigger points, onset and offset
    
    properties
        
        coder
        curSess
        dataType
        dataUnit
        
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
        
        sigsLims
        
        timeSec
        sigsSec
        
        sigsBase
        sigsBaseM
        
        sigsNorm
        sigsNormSec
        
        sigsSecM
        sigsSecMLims
        
        sigsMeanFig
        sigsMeanFigTitle
        
        dataColor1
        dataColor2
        dataColor3
        
        OnsetOffsetAxes
        legendCurves
        legendLabels
        LgdObj
    end
    
    methods
        function obj = dfSectionDataOrg(time, sigs, trigs, fs, dataInfo)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.coder    = dataInfo.coder;
            obj.curSess  = dataInfo.curSess;
            obj.dataType = dataInfo.sigType;
            obj.dataUnit = dataInfo.units;
            
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
            
            % Time vector corresponding to the sectioned signals
            obj.timeSec = linspace(obj.preEveT, obj.posEveT, obj.numSampSec);
            
            obj.dataColor1 = [55,126,184]/255; % Cerulean Blu
            obj.dataColor2 = [77,175,74]/255;  % Leaf Green
            obj.dataColor3 = [231,41,138]/255; % Bright Magenta
            
            obj.legendCurves = [];
            obj.legendLabels = {};
            
%             % Section raw f0 around onset and offset
%             obj.sigsSec  = obj.sectionData(obj.sigs);
            
            % Identify baseline values
%             obj = obj.identifyBaselineValues(obj.sigsSec);
            
            % Convert to cents
%             obj = obj.convertCentsData();
            
            % Section converted f0 around onset and offset
%             obj.sigsNormSec = obj.sectionData(obj.sigsNorm);
            
            % Identify trials to be removed.
            
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
        
        function obj = identifyBounds(obj)
            
            % Full Trial Limits
            uLSigs = max(obj.sigs);
            lLSigs = min(obj.sigs);
            
            uL = round(max(uLSigs)) + 10;
            lL = round(min(lLSigs)) - 10;
            
            obj.sigsLims = [obj.time(1) obj.time(end) lL uL];
            
            % Section Mean Limits
            [lwBoundOn, upBoundOn] = obj.MaxMinBounds(obj.sigsSecM(:,1), obj.sigsSecM(:,2), 10);
            [lwBoundOf, upBoundOf] = obj.MaxMinBounds(obj.sigsSecM(:,3), obj.sigsSecM(:,4), 10);

            lwBoundM = obj.CompareAndChooseBounds(lwBoundOn, lwBoundOf, 'min');
            upBoundM = obj.CompareAndChooseBounds(upBoundOn, upBoundOf, 'max');

            obj.sigsSecMLims = [obj.timeSec(1) obj.timeSec(end) lwBoundM upBoundM];
        end
        
        function [lwBound, upBound] = MaxMinBounds(obj, audio, audioErr, buff)
            % MaxMinBounds takes an audio trace and the error trace for that audio, and
            % identifies what the max and min values of the trace are. It identifies
            % how wide the bounds needs to be fully show the trace with some buff

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
        
        function obj = drawSigsSecM(obj)
            % drawDAQMeanTrialMicf0(res, plotFolder) plots differences in microphone 
            % recordings between perturbed and control trials. 

%             if isempty(varargin)
%                 presFlag = 0;
%             else
%                 presFlag = varargin{1};
%             end
% 
%             curSess          = res.curSess;
%             f0b              = round(res.f0b, 1); % Baseline f0 rounded to 0.1 Hz
%             AudFB            = res.AudFB;
%             numCT            = res.numContTrialsFin;
%             numPT            = res.numPertTrialsFin;
% 
%             time             = res.secTime;
%             meanf0PertOnset  = res.audioMf0MeanPert(:,1);
%             CIf0PertOnset    = res.audioMf0MeanPert(:,2);
%             meanf0PertOffset = res.audioMf0MeanPert(:,3);
%             CIf0PertOffset   = res.audioMf0MeanPert(:,4);
% 
%             meanf0ContOnset  = res.audioMf0MeanCont(:,1);
%             CIf0ContOnset    = res.audioMf0MeanCont(:,2);
%             meanf0ContOffset = res.audioMf0MeanCont(:,3);
%             CIf0ContOffset   = res.audioMf0MeanCont(:,4);
%             limits           = res.limitsAmean;

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
            fontName       = 'Arial';
            titleFontSize  = 14;
            axisLSize      = 14;
            lineThick      = 4;

            % Subplot dimensions and spacing
            obj.OnsetOffsetAxes = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

            % Onset of Perturbation
            axes(obj.OnsetOffsetAxes(1))
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            dHOn = shadedErrorBar(obj.timeSec, obj.sigsSecM(:,1), obj.sigsSecM(:,2), 'lineprops', {'color', obj.dataColor1}, 'transparent', 1);

            set(dHOn.mainLine, 'LineWidth', lineThick)
            xlabel('Time (s)', 'FontName', fontName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', fontName, 'FontSize', axisLSize, 'FontWeight', 'bold')
            title('Onset of Perturbation', 'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off

            set(gca,'FontName', fontName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold')
            hold off
            
            %Offset of Perturbation
            axes(obj.OnsetOffsetAxes(2))   
            plot(zeroLineX, zeroLineY, 'color', zeroLineC, 'LineWidth', lineThick, 'LineStyle', '--')
            hold on
            plot(trigLineX, trigLineY, 'color', trigLineC, 'LineWidth', lineThick)
            dHOf = shadedErrorBar(obj.timeSec, obj.sigsSecM(:,3), obj.sigsSecM(:,4), 'lineprops', {'color', obj.dataColor1}, 'transparent', 1);
              
            obj.legendCurves = cat(2, obj.legendCurves, dHOf.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, ['Line 1: ' num2str(obj.numTrial) ' Trials']);
            
            set(dHOf.mainLine, 'LineWidth', lineThick)
            xlabel('Time (s)', 'FontName', fontName, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
            ylabel([obj.dataType ' (' obj.dataUnit ')'], 'FontName', fontName, 'FontSize', axisLSize, 'FontWeight', 'bold')
            title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
            axis(obj.sigsSecMLims); box off
            hold off
            set(gca,'FontName', fontName,...
                    'FontSize', axisLSize,...
                    'FontWeight','bold',...
                    'YAxisLocation', 'right');

            sup = suptitle({obj.curSess});
            set(sup, 'FontName', fontName,...
                     'FontSize', titleFontSize,...
                     'FontWeight','bold')

%             annoStim = ['SM: ' num2str(statSM) ' cents'];
%             annoResp = ['RM: ' num2str(statRM) ' cents'];
%             annoPerc = ['RP: ' num2str(statRP) ' %'];
% 
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
        
        function obj = appendFigure(obj, sigsSecM, flag)
            
            if flag == 2
                color = obj.dataColor2;
            else
                color = obj.dataColor3;
            end
            
            axes(obj.OnsetOffsetAxes(1))
            hold on
            onH = shadedErrorBar(obj.timeSec, sigsSecM(:,1), sigsSecM(:,2), 'lineprops', {'color', color}, 'transparent', 1);
            set(onH.mainLine, 'LineWidth', 4)
            hold off
            
            axes(obj.OnsetOffsetAxes(2))
            hold on
            ofH = shadedErrorBar(obj.timeSec, sigsSecM(:,3), sigsSecM(:,4), 'lineprops', {'color', color}, 'transparent', 1);
            set(ofH.mainLine, 'LineWidth', 4)
            hold off
            
            obj.legendCurves = cat(2, obj.legendCurves, ofH.mainLine);
            obj.legendLabels = cat(2, obj.legendLabels, ['Line ' num2str(flag) ': ' num2str(obj.numTrial) ' Trials']);
            
            legend(obj.legendCurves, obj.legendLabels);
        end
        
        function saveSigsSecMFig(obj, plotFolder)
            
            saveFileName = fullfile(plotFolder, obj.sigsMeanFigTitle);
            export_fig(saveFileName)
        end
    end
end