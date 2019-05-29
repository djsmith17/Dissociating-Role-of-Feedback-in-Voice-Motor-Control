classdef MeasureSummaryStats
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pAnalysis
        SavResultsDir
        alphaLevel
        SummaryStruct
        SummaryTable
        statSentTable
    end
    
    methods
        function obj = MeasureSummaryStats(dirs, pA, varName, cond, measure, idealLambda)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            numObs    = length(measure);

            % Load necessary identifiers
            obj.pAnalysis     = pA.pAnalysis;
            obj.SavResultsDir = dirs.SavResultsDir;
            obj.alphaLevel    = 0.05/3;
            
            % Unpack the measure into a structure
            str.varName  = varName;
            str.cond     = cond;
            str.measure  = measure;        % Raw Data Values
            
            % Calculate the Descriptive Stats
            str.mean     = round(mean(measure), 2);
            str.median   = round(median(measure), 2);
            str.min      = round(min(measure), 2);
            str.max      = round(max(measure), 2);
            str.SD       = round(std(measure), 2);
            str.SE       = round(str.SD/sqrt(numObs), 2);

            str.isTrans     = 0;       % Default is not transformed
            str.measureT    = measure; % Transformed Data Values (Default is the same)
            str.measureZ    = [];      % Z-Scored Data Values
            str.idealLambda = idealLambda;
            str.usedLambda  = 'N/A';   % Default is not transformed
            str.suffix      = '';      % Default is not transformed

            tbl        = table();
            tbl.mean   = str.mean;
            tbl.min    = str.min;
            tbl.median = str.median;
            tbl.max    = str.max;
            tbl.SD     = str.SD;
            tbl.SE     = str.SE;
            tbl.Properties.RowNames = {cond};
            
            obj.SummaryStruct = str;
            obj.SummaryTable  = tbl;
        end
        
        function obj = testNormality(obj)
            % Skew and Kurtosis
            obj.SummaryStruct.measureSkew     = round(skewness(obj.SummaryStruct.measureT), 4);
            obj.SummaryStruct.measureKurtosis = round(kurtosis(obj.SummaryStruct.measureT), 2);

            % Z-Score and Shapiro-Wilk Test
            obj.SummaryStruct.measureZ  = zscore(obj.SummaryStruct.measureT);
            [swH, swPValue, swTest]     = swtest(obj.SummaryStruct.measureZ);

            % Add to the Summmary Data Structure
            obj.SummaryStruct.swH      = double(swH);
            obj.SummaryStruct.swPValue = round(swPValue, 3);
            obj.SummaryStruct.swTest   = round(swTest, 3);

            % Populate Summary Stat Table
            obj.SummaryTable.Skew     = obj.SummaryStruct.measureSkew;
            obj.SummaryTable.Kurtosis = obj.SummaryStruct.measureKurtosis;
            obj.SummaryTable.swH      = obj.SummaryStruct.swH;
            obj.SummaryTable.swPValue = obj.SummaryStruct.swPValue;
            obj.SummaryTable.swTest   = obj.SummaryStruct.swTest;
        end
        
        function obj = performSimpleBoxCoxTrans(obj)
            
            [obj.SummaryStruct.measureT, l] = boxcox(obj.SummaryStruct.measure + 1 - obj.SummaryStruct.min);
            obj.SummaryStruct.isTrans    = 1;
            obj.SummaryStruct.suffix     = 'Trans';
            obj.SummaryStruct.usedLambda = num2str(round(l,2));
            
        end
        
        function obj = performTTest(obj)
            % Perform a One-Sample T-Test on the difference between the measures
            
            [~, P, ~, STATS] = ttest(obj.SummaryStruct.measureT);
            Pstr   = sprintf('%0.6f', P);
            
            if P < obj.alphaLevel
                isSig = 1;
                sigNote = '';
            else
                isSig = 0;
                sigNote = ' not';
            end
            
            statSentence = sprintf('The variable %s was%s significantly different between the two conditions t(%d) = %0.2f, p = %0.6f\n',...
                                   obj.SummaryStruct.varName,...
                                   sigNote,...
                                   STATS.df,...
                                   STATS.tstat,...
                                   P);
            
            obj.SummaryStruct.ttestStat   = STATS.tstat;
            obj.SummaryStruct.ttestP      = P;
            obj.SummaryStruct.ttestPstr   = Pstr;
            obj.SummaryStruct.isSig       = isSig;
            obj.SummaryStruct.statSentence = statSentence;
            obj.statSentTable = table({statSentence}, 'VariableNames', {'StatSentence'});
            
        end
        
        function drawHistoBoxCombo(obj)

            measure    = obj.SummaryStruct.measureT;
            varName    = obj.SummaryStruct.varName;
            suffix     = obj.SummaryStruct.suffix;
            usedLambda = obj.SummaryStruct.usedLambda;
            swH        = num2str(obj.SummaryStruct.swH); 
            swP        = num2str(obj.SummaryStruct.swPValue); 
            swW        = num2str(obj.SummaryStruct.swTest);
            
            measureM   = num2str(obj.SummaryStruct.mean);
            measureSD  = num2str(obj.SummaryStruct.SD);

            lambda = '\lambda';
            mu = '\mu';
            sigma  = '\sigma'; 

            diffBox = figure('Color', [1 1 1]);
            plotpos = [30 0]; plotdim = [800 300];
            set(diffBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

            subplot(1,2,1); histogram(measure, 10); box off
            title(['H=' swH ', p=' swP ', W=' swW])

            subplot(1,2,2); boxplot(measure); box off
            suptitle(varName)

            annotation('textbox',[0.8 0.88 0.45 0.1],...
                       'string', {[lambda ' = ' usedLambda]},...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',14,...
                       'FontName','Arial');
                   
            annotation('textbox',[0.80 0.48 0.45 0.1],...
                       'string', {[mu ' = ' measureM 'psi'],...
                                  [sigma ' = ' measureSD 'psi']},...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',14,...
                       'FontName','Arial');

            BoxPlotFigureFile = fullfile(obj.SavResultsDir, [obj.pAnalysis varName suffix 'BoxPlotCombo.jpg']);
            export_fig(BoxPlotFigureFile)
        end
    end
end