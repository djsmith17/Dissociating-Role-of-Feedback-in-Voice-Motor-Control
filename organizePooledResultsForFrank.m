function organizePooledResultsForFrank(dirs, allSubjRes)

time          = allSubjRes.secTime';
pertResponses = allSubjRes.audioMf0SecPert{1};

onsets  = pertResponses(:, :, 1);
offsets = pertResponses(:, :, 2);

onsetsRatio  = exp(onsets/1200);
offsetsRatio = exp(offsets/1200);

onsetSet  = [time, onsetsRatio];
offsetSet = [time, offsetsRatio];

dirs.alignedMicResponsesFile = fullfile(dirs.SavResultsDir, 'Pitch-Shift Reflex Aligned Mic Recordings.xlsx');
xlswrite(dirs.alignedMicResponsesFile, onsetSet, 'Perturbation Onset');
xlswrite(dirs.alignedMicResponsesFile, offsetSet, 'Perturbation Offset');
end