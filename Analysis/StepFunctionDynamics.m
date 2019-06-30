classdef StepFunctionDynamics
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time  
        sensor
        fs
        
        pertIdx
        pertTime
        
        curSess
        
        numTrial

        % Per Trial sensor index and time value of rising edge start and falling
        % edge end (expecting ~step function)
        TrigIdx  
        TrigTime 

        % Per Trial lag time of trigger to start of rising edge
        lagTimes  
        lagTimeM  
        lagTimeSE

        % Per Trial rise time of rising edge (step function)
        riseTimes
        riseTimeM
        riseTimeSE

        % Per trial value at Onset/Offset (pressure, voltage, etc)
        OnOffVal
        OnOffValM
        OnOffValSE

        % Per trial loss in value (pressure, voltage, etc)
        pTrialLoss
        pTrialLossM
        pTrialLossSE
        
        % Sensor recording aligned at perturbation trigger
        timeAl
        sensorAl

        % Sensor recording sectioned around onset/offset
%         timeSec
%         sensorSec
%         sensorSecM

        sensorLim
        sensorLimAl
    end
    
    methods
        function obj = StepFunctionDynamics(time, sensor, fs, pertIdx, pertTime, curSess)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

            % Analyzing the dynamics of the sensor during onset and offset of the
            % stimulus signal.
            
            obj.time   = time;
            obj.sensor = sensor;
            obj.fs     = fs;
            
            obj.pertIdx  = pertIdx;
            obj.pertTime = pertTime;
            
            obj.curSess = curSess;

            [~, obj.numTrial] = size(obj.sensor);

            for ii = 1:obj.numTrial
                trial = obj.sensor(:,ii);

                % 1st Derivative (1D) of the recording
                fDiff = [0; diff(trial)];
                fDiff = smooth(5*fDiff, 20);

                % Thresholds for detecting edges of (expected) step function
                threshUp = 0.2*max(fDiff);
                threshDn = 0.3*min(fDiff);

                ups = find(fDiff > threshUp); % Parts of 1D that could include Increasing edge
                dns = find(fDiff < threshDn); % Parts of 1D that could include decreasing edge

                StRiseIdx = ups(1);           % Assume first idx of increasing edges is (St)art of rise 
                StFallIdx = dns(1);           % Assume first idx of decreasing edges is (St)art of fall

                RisingEdgeRange = StRiseIdx:(StRiseIdx + 0.3*fs); % Range Following the StRiseIdx
                [~, idxAtMax] = max(trial(RisingEdgeRange));      % Max value of recording in that range (Assume low->high)
                SpRiseIdx = StRiseIdx + idxAtMax-1;               % Idx of Rise (S)to(p)

                FallingEdgeRange = StFallIdx:(StFallIdx + 0.3*fs); % Range Following the StFallIdx
                [~, idxAtMin] = min(trial(FallingEdgeRange));      % Min value of recording in that range (Assume high->low)
                SpFallIdx = StFallIdx + idxAtMin-1;                % Idx of Fall (S)to(p)

                % Convert Indices to Times
                StRiseTime = round(time(StRiseIdx), 3);
                SpRiseTime = round(time(SpRiseIdx), 3);
                StFallTime = round(time(StFallIdx), 3);
                SpFallTime = round(time(SpFallIdx), 3);

                lagTimeRise = StRiseTime - pertTime(ii, 1); 
                lagTimeFall = StFallTime - pertTime(ii, 2);

                riseTime = SpRiseTime - StRiseTime;
                fallTime = SpFallTime - StFallTime;

                % What was the val of the sensor at the end of Rising Edge (Start Plateau)
                RiseVal = trial(SpRiseIdx);

                % What was the val of the sensor at the start of the Falling Edge (End Plateau)
                FallVal = trial(StFallIdx);

                obj.TrigIdx  = cat(1, obj.TrigIdx, [StRiseIdx StFallIdx]);
                obj.TrigTime = cat(1, obj.TrigTime, [StRiseTime StFallTime]);

                obj.lagTimes  = cat(1, obj.lagTimes, [lagTimeRise lagTimeFall]);
                obj.riseTimes = cat(1, obj.riseTimes, [riseTime fallTime]);

                obj.OnOffVal = cat(1, obj.OnOffVal, [RiseVal, FallVal]);
            end
            
            %Difference in val at between subsequent trials at SpRiseIdx
            obj.pTrialLoss = diff(obj.OnOffVal(:,1));

            %Mean and SE lag time
            obj.lagTimeM  = mean(obj.lagTimes);
            obj.lagTimeSE = (std(obj.lagTimes))/sqrt(obj.numTrial);

            %Mean and SE rise time
            obj.riseTimeM  = mean(obj.riseTimes);
            obj.riseTimeSE = (std(obj.riseTimes))/sqrt(obj.numTrial);

            %Mean and SE Val at SpRiseIdx
            obj.OnOffValM  = mean(obj.OnOffVal);
            obj.OnOffValSE = (std(obj.OnOffVal))/sqrt(obj.numTrial);

            %Mean and SE trial-trial difference of val at SpRiseIdx
            obj.pTrialLossM  = mean(obj.pTrialLoss);
            obj.pTrialLossSE = (std(obj.pTrialLoss))/sqrt(obj.numTrial-1);

            % Round and Covert to standard units
            obj.lagTimeM  = round(1000*obj.lagTimeM, 3); % Round and Convert s -> ms
            obj.lagTimeSE = round(1000*obj.lagTimeSE, 3); % Round and Convert s -> ms

            obj.riseTimeM  = round(1000*obj.riseTimeM, 3); % Round and Convert s -> ms
            obj.riseTimeSE = round(1000*obj.riseTimeSE, 3); % Round and Convert s -> ms

            obj.OnOffValM  = round(obj.OnOffValM, 3); % Round
            obj.OnOffValSE = round(obj.OnOffValSE, 3); % Round

            obj.pTrialLossM  = round(obj.pTrialLossM, 3); % Round
            obj.pTrialLossSE = round(obj.pTrialLossSE, 3); % Round

            %%%%%%%%%%%%%
            [obj.timeAl, obj.sensorAl] = obj.alignSensorData(obj.sensor, obj.fs, pertIdx);

%             [obj.timeSec, obj.sensorSec] = sectionData(sensor, fs, pertIdx);
%             obj.sensorSecM               = meanSensorData(obj.sensorSec);

            maxPres = max(max(sensor));
            minPres = min(min(sensor));
            uLP = maxPres + 0.2;
            lLP = minPres - 0.2;

            obj.sensorLim   = [0 4 lLP uLP];

            %Aligned Pressure Data
            obj.sensorLimAl = [-1 2.5 lLP uLP];
        end
        
        function [timeAl, sensorAl] = alignSensorData(obj, sensor, fs, idx)
        % [timeAl, sensorAl]  = alignSensorData(sensor, fs, idx) sections 
        % multi-trial sensor data about individual trial trigger points. 
        % Each sectioned trial is of equal length, and includes equal lengths of 
        % data on either side of the trigger point. 
        % The sectioned trials are then concatenated into a matrix, which are
        % aligned the trigger point of each trial. 
        %
        % sensor: recorded sensor data (numSamp x numTrial)
        % fs:     sampling rate of sensor data
        % idx:    Onset and Offset trigger POINTS (numTrial, 2)
        %
        % timeAl:   vector of time points corresponding to the sectioned data
        % sensorAl: sectioned and aligned sensor data 

        preEve = 0.5; % time preEvent Seconds 
        posEve = 2.0; % time posEvent Seconds

        % At the moment only aligning by Onset. 
        % This could eventually become an input to toggle between Onset/Offset
        OnsetTrigs = idx(:, 1); 

        sensorAl = [];
        for ii = 1:obj.numTrial
            St = OnsetTrigs(ii) - fs*preEve;  % Points preEvent (trigger)
            Sp = OnsetTrigs(ii) + fs*posEve;  % Points posEvent (trigger)

            sensorSec = sensor(St:Sp, ii); % Grab St:Sp around the trigger point for this trial

            sensorAl = cat(2, sensorAl, sensorSec);
        end
        per     = 1/fs;

        % Time Vector of the sectioned data
        timeAl = (-preEve:per:posEve)';
        end
        
        function drawSingleTrial(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            figure
            plot(obj.time, trial)
            hold on
            plot([StRiseTime StRiseTime], [-5 10], 'g')
            hold on
            plot([SpRiseTime SpRiseTime], [-5 10], 'g--')
            hold on
            plot([StFallTime StFallTime], [-5 10], 'r')
            hold on
            plot([SpFallTime SpFallTime], [-5 10], 'r--')
            axis([0 4 -0.5 5])
        end
        
        function drawAllTrial(obj)
            
            plotpos = [500 300]; plotdim = [2000 1000];
            allTrialFig = figure('Color', [1 1 1]);
            set(allTrialFig, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

            onsetColor     = 'g';
            offsetColor    = 'r';
            pertBoxC     = [0.8 0.8 0.8];
            fontN        = 'Arial';
            legAnnoFSize = 12;
            titleFSize   = 14;
            axisLSize    = 14;
            lineThick    = 4;

            numColums = 5;
            numrows   = ceil(obj.numTrial/numColums);

            ha = tight_subplot(numrows, numColums, [0.08 0.03],[0.05 0.05],[0.05 0.03]);

            for ii = 1:obj.numTrial      
                axes(ha(ii))
                
                pertAx  = [obj.pertTime(ii,1), obj.pertTime(ii,2)];
                pertAy  = [500 500];
                
                area(pertAx, pertAy, -600, 'FaceColor', pertBoxC, 'EdgeColor', pertBoxC)
                hold on
                
                plot(obj.time, obj.sensor(:,ii), 'k', 'LineWidth', 1.5)
                hold on
                plot([obj.TrigTime(ii,1) obj.TrigTime(ii,1)], [-5 10], 'Color', onsetColor)
                hold on
                plot([obj.TrigTime(ii,1)+obj.riseTimes(ii,1) obj.TrigTime(ii,1)+obj.riseTimes(ii,1)], [-5 10], 'Color', onsetColor, 'LineStyle', '--')
                hold on
                plot([obj.TrigTime(ii,2) obj.TrigTime(ii,2)], [-5 10], 'Color', offsetColor)
                hold on
                plot([obj.TrigTime(ii,2)+obj.riseTimes(ii,2) obj.TrigTime(ii,2)+obj.riseTimes(ii,2)], [-5 10], 'Color', offsetColor, 'LineStyle', '--')

                if ii == 1
                    xlabel('Time (s)')
                    ylabel('Pressure (psi)')
                end
                axis(obj.sensorLim);
                box off
                set(gca,'FontSize', 14,...
                        'FontWeight','bold')
                    

%                 title({['Trial ' num2str(trialNums(ii))], [num2str(pertTrig(ii,1)) '  ' num2str(pertTrig(ii,2))]}, 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
%                 axis(limits); box off
% 
%                 set(gca, 'FontName', fontN,...
%                          'FontSize', axisLSize,...
%                          'FontWeight','bold')
            end
            suptitle(obj.curSess)
            
            dirPath = 'E:\Desktop\Pressure Diagnostics';
            filePath = fullfile(dirPath, [obj.curSess 'PressureDiag.jpg']);
            export_fig(filePath)
            
        end
        
        function drawAlignedSensor(obj)
        %Plots multiple trials on top of each other. Currently only plotting one 
        %sensor. Assumes the trials have been aligned.


        LagNote = ['Mean Onset Lag: ' num2str(obj.lagTimeM(1)) ' ms ± ' num2str(obj.lagTimeSE(1)) ' ms'];

        RiseNote = ['Mean Rise Time: ' num2str(obj.riseTimeM(1)) ' ms ± ' num2str(obj.riseTimeSE(1)) ' ms'];

        PlatStNote = ['Mean Plataeu (Start) Val: ' num2str(obj.OnOffValM(1)) ' psi ± ' num2str(obj.OnOffValSE(1)) ' psi'];

        PresLossNote = ['Mean Pressure Loss: ' num2str(obj.pTrialLossM(1)) ' psi ± ' num2str(obj.pTrialLossM(1)) ' psi'];

        plotpos = [500 300];
        plotdim = [850 600];
        trialColors = distinguishable_colors(obj.numTrial);

        CombinedSensor = figure('Color', [1 1 1]);
        set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

        plot([0 0], [-1 20], 'k-', 'LineWidth', 2)

        for ii = 1:obj.numTrial
            hold on
            h(ii) = plot(obj.timeAl, obj.sensorAl(:,ii), 'LineWidth', 2, 'Color', trialColors(ii,:));
        end

        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
        ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
        title('Mean Pressure Sensor Measurements aligned at Perturbation Onset', 'FontSize', 12, 'FontWeight', 'bold')
        axis(obj.sensorLimAl);
        box off

        set(gca,'FontSize', 12,...
                'FontWeight', 'bold')

        lgdNames = cell(obj.numTrial, 1);
        for i = 1:obj.numTrial; lgdNames{i} = ['Trial ' num2str(i)]; end

        pltlgd = legend(h, lgdNames);
        set(pltlgd, 'box', 'off',...
                    'location', 'NorthWest'); 

        t = annotation('textbox',[0.64 0.83 0.9 0.1],...
                       'string', {LagNote;...
                                  RiseNote;...
                                  PlatStNote;...
                                  PresLossNote},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',9,...
                        'FontName','Arial');
        end
    end
end

