function organizePooledResultsForFrank(dirs, allSubjRes)

time          = allSubjRes.secTime';
pertResponses = allSubjRes.audioMf0SecPert{1};

onsets  = pertResponses(:, :, 1);
offsets = pertResponses(:, :, 2);

onsetsRatio  = exp(onsets/1200);
offsetsRatio = exp(offsets/1200);

onsetSet  = [time, onsetsRatio];
offsetSet = [time, offsetsRatio];

dirs.alignedMicResponsesFile_Onset  = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Onset.csv');
dirs.alignedMicResponsesFile_Offset = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings_Offset.csv');
csvwrite(dirs.alignedMicResponsesFile_Onset, onsetSet);
csvwrite(dirs.alignedMicResponsesFile_Offset, offsetSet);
end