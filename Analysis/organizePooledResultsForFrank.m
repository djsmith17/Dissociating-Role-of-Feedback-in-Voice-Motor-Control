function organizePooledResultsForFrank(dirs, allSubjRes)

% drag these out from the pooled result structures
time          = round(allSubjRes.secTime', 4);
pertResponses = allSubjRes.audioMf0SecPert{1};

downSamp = find(rem(time, 0.01) == 0);

timeDN          = time(downSamp);
pertResponsesDN = pertResponses(downSamp, :,:);

% separate the onset and offset responses
onsets        = pertResponsesDN(:, :, 1);
offsets       = pertResponsesDN(:, :, 2);

% recreate the stimulus
tStep = '10ms';
t0  = 0;
tRD = 0.110; % 110ms down
tRU = 0.150; % 150ms up

zeroPt   = find(timeDN == t0);
rampDnPt = find(timeDN == tRD);
rampUpPt = find(timeDN == tRU);

rampDnLen = length(zeroPt:rampDnPt);
rampUpLen = length(zeroPt:rampUpPt);

finalpert = 2^(-100/1200);

stim = ones(length(timeDN), 1);

stimOnset = stim;
stimOnset(rampDnPt:end) = finalpert;
stimOnset(zeroPt:rampDnPt) = linspace(1, finalpert, rampDnLen);
stimOnset = stimOnset -1;

stimOffset = stim;
stimOffset(1:zeroPt) = finalpert;
stimOffset(zeroPt:rampUpPt) = linspace(finalpert, 1, rampUpLen);
stimOffset = stimOffset -1;

onsetSet  = [stimOnset, onsets];
offsetSet = [stimOffset, offsets];

dirs.alignedMicResponsesFile_Onset  = fullfile(dirs.SavResultsDir, ['Pitch-Shift Reflex Aligned Mic Recordings_Onset' tStep '.csv']);
dirs.alignedMicResponsesFile_Offset = fullfile(dirs.SavResultsDir, ['Pitch-Shift Reflex Aligned Mic Recordings_Offset' tStep '.csv']);
csvwrite(dirs.alignedMicResponsesFile_Onset, onsetSet);
csvwrite(dirs.alignedMicResponsesFile_Offset, offsetSet);
end