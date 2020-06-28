function [IntraCoderSentence, IntraCoderDiffSentence] = dfIdentifyIntraCoderReliabilityScore()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

eAn.allParti = {'DRF5', 'DRF14', 'DRF19'};
eAn.numParti = length(eAn.allParti);
eAn.runs     = 'SFL1';
eAn.trials   = {'3' '10' '7'};
eAn.trialIdx = [3 9 5];
eAn.coders   = {'RF', 'RF2'};

varNames = {'CurSessTrial', 'Frame', 'Line', 'EuclidianDist1', 'EuclidianDist2', 'Diff'};
varTypes = {'string', 'double' 'double', 'double', 'double', 'double'};
numVar = length(varNames);
reliTable = table('Size', [1 numVar], 'VariableTypes', varTypes, 'VariableNames', varNames);

tableRow = 1;
for ii = 1:eAn.numParti
    participant = eAn.allParti{ii};
    run         = eAn.runs;
    trialName   = ['Trial ' eAn.trials{ii}];
    trialIdx    = eAn.trialIdx(ii);
    
    coder1      = eAn.coders{1};
    coder2      = eAn.coders{2};
    
    dirs.ResultsParti     = fullfile(dirs.Results, participant, run);
    dirs.ResultsCodedVid1 = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasures' coder1 '.mat']);
    dirs.ResultsCodedVid2 = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasures' coder2 '.mat']);
    
    load(dirs.ResultsCodedVid1); % returns CodedEndoFrameDataSet
    CodedSet1 = CodedEndoFrameDataSet;
    load(dirs.ResultsCodedVid2); % returns CodedEndoFrameDataSet
    CodedSet2 = CodedEndoFrameDataSet;
    
    curTable1 = CodedSet1{trialIdx};
    curTable2 = CodedSet2{trialIdx};
    [numFrames, ~] = size(curTable1);
    numLines = 3;
    
    for jj = 1:numFrames
        for kk = 1:numLines
            reliTable.CurSessTrial(tableRow) = [participant run '_' trialName];
            reliTable.Frame(tableRow)        = jj;
            reliTable.Line(tableRow)         = kk;
            switch kk
                case 1
                    reliTable.EuclidianDist1(tableRow) = curTable1.Dist(jj);
                    reliTable.EuclidianDist2(tableRow) = curTable2.Dist(jj);
                case 2
                    reliTable.EuclidianDist1(tableRow) = curTable1.Dist2(jj);
                    reliTable.EuclidianDist2(tableRow) = curTable2.Dist2(jj);
                case 3
                    reliTable.EuclidianDist1(tableRow) = curTable1.Dist3(jj);
                    reliTable.EuclidianDist2(tableRow) = curTable2.Dist3(jj);
            end
            
            reliTable.Diff(tableRow) = reliTable.EuclidianDist2(tableRow) - reliTable.EuclidianDist1(tableRow);
            
            tableRow = tableRow + 1;
        end   
    end
end

IntraCoderResposne  = [reliTable.EuclidianDist1 reliTable.EuclidianDist2];
[corrR, corrP] = corrcoef(IntraCoderResposne);

IntraCoderSentence = sprintf('Strong correlation (r = %0.3f, p = %0.3f), between Coding sessions 1 and 2, (n = %d).', corrR(2,1), corrP(2,1), length(IntraCoderResposne));

diffMean = mean(reliTable.Diff);
diffSTD  = std(reliTable.Diff);

IntraCoderDiffSentence = sprintf('The mean difference in measured euclidian distance between the two coding sessions was %0.2f pixels (SD = %0.2f pixels).', diffMean, diffSTD);
end