function organizePooledResultsForFrank(dirs, allSubjRes)

% drag these out from the pooled result structures
time          = round(allSubjRes.secTime', 4);
pertResponses = allSubjRes.audioMf0SecPert{1};

% separate the onset and offset responses
onsets        = pertResponses(:, :, 1);
offsets       = pertResponses(:, :, 2);

% Create a pert
t0  = 0;
tRD = 0.110; % 110ms down
tRU = 0.150; % 150ms up
finalpert = 2^(-100/1200); % 100 cents converted to Hz 

zeroPt   = find(time == t0);
rampDnPt = find(time == tRD);
rampUpPt = find(time == tRU);

rampDnLen = length(zeroPt:rampDnPt);
rampUpLen = length(zeroPt:rampUpPt);

stim = ones(length(time), 1);

stimOnset = stim;
stimOnset(rampDnPt:end) = finalpert;
stimOnset(zeroPt:rampDnPt) = linspace(1, finalpert, rampDnLen);
stimOnset = stimOnset -1;

stimOffset = stim;
stimOffset(1:zeroPt) = finalpert;
stimOffset(zeroPt:rampUpPt) = linspace(finalpert, 1, rampUpLen);
stimOffset = stimOffset -1;

%Align signals
onsetSet  = [stimOnset, onsets];
offsetSet = [stimOffset, offsets];

% Perform the downsample
tStep = '10ms';
downSamp = find(rem(time, 0.01) == 0);
onsetSetDN  = onsetSet(downSamp, :);
offsetSetDN = offsetSet(downSamp, :);

dirs.alignedMicResponsesFile_Onset  = fullfile(dirs.SavResultsDir, ['Pitch-Shift Reflex Aligned Mic Recordings_Onset' tStep 'Hz.csv']);
dirs.alignedMicResponsesFile_Offset = fullfile(dirs.SavResultsDir, ['Pitch-Shift Reflex Aligned Mic Recordings_Offset' tStep 'Hz.csv']);
csvwrite(dirs.alignedMicResponsesFile_Onset, onsetSetDN);
csvwrite(dirs.alignedMicResponsesFile_Offset, offsetSetDN);
end