function organizePooledResultsForFrank(dirs, allSubjRes)

time = allSubjRes.secTime';
zeroPt = 501;
rampDnPt = 611;
rampUpPt = 651;

rampDnLen = length(zeroPt:rampDnPt);
rampUpLen = length(zeroPt:rampUpPt);

finalpert = 2^(-100/1200);

stim = ones(length(time), 1);

stimOnset = stim;
stimOnset(rampDnPt:end) = finalpert;
stimOnset(zeroPt:rampDnPt) = linspace(1, finalpert, rampDnLen);
stimOnset = stimOnset -1;

stimOffset = stim;
stimOffset(1:zeroPt) = finalpert;
stimOffset(zeroPt:rampUpPt) = linspace(finalpert, 1, rampUpLen);
stimOffset = stimOffset -1;

pertResponses = allSubjRes.audioMf0SecPert{1};

onsets  = pertResponses(:, :, 1);
offsets = pertResponses(:, :, 2);

onsetSet  = [stimOnset, onsets];
offsetSet = [stimOffset, offsets];

dirs.alignedMicResponsesFile_Onset  = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Onset.csv');
dirs.alignedMicResponsesFile_Offset = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Offset.csv');
csvwrite(dirs.alignedMicResponsesFile_Onset, onsetSet);
csvwrite(dirs.alignedMicResponsesFile_Offset, offsetSet);
end