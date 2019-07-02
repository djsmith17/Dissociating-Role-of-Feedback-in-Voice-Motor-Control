classdef MeasureSummaryStats
    % MeasureSummaryStats(dirs, pA, measureVar, measure, idealLambda) 
    % is a class for organizing the resultant stats of a measure. The input
    % dirs is the current directory you are working in and
    % specifically is looking for the the results folder you are
    % writing to. pA is a structure describing the pooled analysis
    % conditions we are testing over. varName is the name of the
    % measure we are testing. Cond is the name of the condition
    % that this variable was tested on. Measure is a vector of the
    % actual measured values. idealLambda is a value of lambda such
    % that when measure is transformed using a box-cox transform 
    % with this lambda, the transformed measure is as normal is
    % possible.
    
    properties
        pAnalysis
        SavResultsDir
        alphaLevel
        SummaryStruct
        SummaryTable
        statSentTable
    end
    
    methods
        function obj = MeasureSummaryStats(dirs, pA, measVar, measure, idealLambda)
            % This class is called from the following Stats functions
            %
            % -StatsOrg_MaskingNoiseStudy
            % -StatsOrg_DRF_Som
            % -StatsOrg_DRF_Aud
            % -StatsOrg_DRF_Som_Aud
            % -drawExpPressureDist

            % Load necessary identifiers
            obj.pAnalysis     = pA.pAnalysis;
            obj.SavResultsDir = dirs.SavResultsDir;
            obj.alphaLevel    = 0.05/3;
            
            % Unpack the measure into a structure
            str.varName  = measVar.varName;
            str.cond     = measVar.condition;
            str.units    = measVar.units;
            str.measure  = measure;        % Raw Data Values
            
            % Calculate the Descriptive Stats
            str.numObvs  = length(measure);
            str.mean     = round(mean(measure), 2);
            str.median   = round(median(measure), 2);
            str.min      = round(min(measure), 2);
            str.max      = round(max(measure), 2);
            str.SD       = round(std(measure), 2);
            str.SE       = round(str.SD/sqrt(str.numObvs), 2);

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
            tbl.Properties.RowNames = {str.cond};
            
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
        
        function obj = performSpecificBoxCoxTrans(obj, l)
            obj.SummaryStruct.measureT   = boxcox(l, obj.SummaryStruct.measure + 1 - obj.SummaryStruct.min);
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
        
        function drawHistogram(obj)
            
            measure    = obj.SummaryStruct.measureT;
            varName    = obj.SummaryStruct.varName;
            units      = obj.SummaryStruct.units;
            suffix     = obj.SummaryStruct.suffix;
            
            numObvs    = num2str(obj.SummaryStruct.numObvs);
            measureM   = num2str(obj.SummaryStruct.mean);
            measureSD  = num2str(obj.SummaryStruct.SD);
            
            mu = '\mu';
            sigma  = '\sigma';
            
            diffBox = figure('Color', [1 1 1]);
            plotpos = [30 40]; plotdim = [400 300];
            set(diffBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

            histogram(measure, 10); box off        
            xlabel({[varName ' (' units ')'], ['(n = ' numObvs ')']})
            
            set(gca, 'FontSize', 12,...
                     'FontWeight', 'bold')

            annotation('textbox',[0.62 0.75 0.45 0.1],...
                       'string', {[mu ' = ' measureM units],...
                                  [sigma ' = ' measureSD units]},...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',12,...
                       'FontName','Arial');

            dirs.DistributionFigureFile = fullfile(obj.SavResultsDir, [obj.pAnalysis varName suffix 'Histogram.jpg']);
            export_fig(dirs.DistributionFigureFile)
        end
        
        function drawHistoBoxCombo(obj)

            measure    = obj.SummaryStruct.measureT;
            varName    = obj.SummaryStruct.varName;
            units      = obj.SummaryStruct.units;
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
            plotpos = [30 40]; plotdim = [800 300];
            set(diffBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

            subplot(1,2,1); histogram(measure, 10); box off
            title(['H=' swH ', p=' swP ', W=' swW])
            
            set(gca, 'FontSize', 12,...
                     'FontWeight', 'bold')

            subplot(1,2,2); boxplot(measure); box off
            suptitle(varName)
            
            set(gca, 'FontSize', 12,...
                     'FontWeight', 'bold')

            % Annotation for Box-Cox Lambda
            annotation('textbox',[0.8 0.88 0.45 0.1],...
                       'string', {[lambda ' = ' usedLambda]},...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',14,...
                       'FontName','Arial');
            
            % Annotation for Mean and SD values    
            annotation('textbox',[0.76 0.64 0.45 0.1],...
                       'string', {[mu ' = ' measureM units],...
                                  [sigma ' = ' measureSD units]},...
                       'LineStyle','none',...
                       'FontWeight','bold',...
                       'FontSize',14,...
                       'FontName','Arial');

            BoxPlotFigureFile = fullfile(obj.SavResultsDir, [obj.pAnalysis varName suffix 'BoxPlotCombo.jpg']);
            export_fig(BoxPlotFigureFile)
        end
    end
end