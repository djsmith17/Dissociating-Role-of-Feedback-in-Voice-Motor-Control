function organizePooledResultsForFrank(dirs, allSubjRes)

time = allSubjRes.secTime';
zeroPt = 501;
rampDnPt = 611;
rampUpPt = 651;

rampDnLen = length(zeroPt:rampDnPt);
rampUpLen = length(zeroPt:rampUpPt);

finalpert = 2^(-100/1200);

set = ones(length(time), 1);

setOnset = set;
setOnset(rampDnPt:end) = finalpert;
setOnset(zeroPt:rampDnPt) = linspace(1, finalpert, rampDnLen);

setOffset = set;
setOffset(1:zeroPt) = finalpert;
setOffset(zeroPt:rampUpPt) = linspace(finalpert, 1, rampUpLen);


pertResponses = allSubjRes.audioMf0SecPert{1};

onsets  = pertResponses(:, :, 1);
offsets = pertResponses(:, :, 2);

onsetsRatio  = 2.^(onsets/1200);
offsetsRatio = 2.^(offsets/1200);

onsetSet  = [setOnset, onsetsRatio];
offsetSet = [setOffset, offsetsRatio];

dirs.alignedMicResponsesFile_Onset  = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Onset.csv');
dirs.alignedMicResponsesFile_Offset = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Offset.csv');
csvwrite(dirs.alignedMicResponsesFile_Onset, onsetSet);
csvwrite(dirs.alignedMicResponsesFile_Offset, offsetSet);
end