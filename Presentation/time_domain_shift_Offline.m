function time_domain_shift_Offline(isOnline, varargin)
%% Offline and online demo for time-domain pitch shifting.
% Usage examples:
% - For offline demo:
%   time_domain_shift_demo('female', 0, '--wav', 'C:\speechres\samples\dante\Pilot22SF2Trial5_micIn.wav', '--play');
%   time_domain_shift_demo('male', 0, '--wav', 'C:\speechres\samples\dante\Pilot26SF3Trial4_micIn.wav', '--play');
%  
% - For online demo:
%   time_domain_shift_demo('female', 1);
%   time_domain_shift_demo('male', 1);

%% Input parameter sanity check and processing.
close all
subject    = 'Pilot24';
run        = 'SF4';
gender     = 'male';  % "male" or "female"
BaseRun    = 'SF3';

expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Auditory Perturbation_Perceptual';
expParam.subject      = subject;
expParam.run          = run;
expParam.gender       = gender;

expParam.baseRun      = BaseRun;
expParam.baseFile     = [expParam.subject expParam.baseRun '.wav'];

toDraw = 1;
toPlay = 0;

%Set our dirs based on the project
dirs = dfDirs(expParam.project);

dirs.SavData     = fullfile(dirs.helpers, 'speechres\samples');
dirs.SavBaseFile = fullfile(dirs.SavData, 'baseline', expParam.baseFile);
dirs.RecData     = fullfile(dirs.SavData, 'data');
dirs.RecWaveDir  = dirs.RecData; 

if ~exist(dirs.SavBaseFile, 'file')
    fprintf('ERROR: No voice file at %s!\n', dirs.SavBaseFile)
    return
end

if isempty(fsic({'female', 'male'}, gender))
    error('Invalid value in pitchGender: %s', gender);
end

expParam.sRate        = 48000;
expParam.frameLen     = 192;
expParam.downFact     = 3;
expParam.sRateAnal    = expParam.sRate/expParam.downFact;
expParam.frameLenDown = expParam.frameLen/expParam.downFact;

%% Get parameters for time-domain shifting.
params = getAudapterDefaultParams(gender);

% Activate time-domain pitch shifting.
params.bTimeDomainShift = 1;

% Set the a priori pitch range.
% Adjust if necessary. You may even consider set this on a per-subject
% basis! The more accurate the a prior estimate, the lower the likelihood
% of incorrect real-time pitch tracking leading to occasional artifacts in
% the pitch shifted feedback.
if isequal(lower(gender), 'female')
    params.pitchLowerBoundHz = 150;
    params.pitchUpperBoundHz = 300;
elseif isequal(lower(gender), 'male')
    params.pitchLowerBoundHz = 80;
    params.pitchUpperBoundHz = 160;
end

% Use a long-enough frame length (and hence spectral/cepstral) window to
% ensure accurate real-time pitch tracking.
params.frameLen = 64;
params.nDelay = 7;
% Time-domain shift requires bCepsLift = 1.
params.bCepsLift = 1;

% The following line specifies the time-domain pitch perturbation schedule.
% It is a 3x2 matrix. Each row corresponds to an "anchor point" in time.
% The 1st element of the row is the time relative to the crossing of the
% intensity threshold (i.e., params.rmsThresh below).
% The 2nd element is the pitch shift ratio. A pitch shift ratio value of
% 1.0 means unshifted. A value of 2 means an octave of up-shift. So 1.0595
% is approximately an upshift of 1 semitone (100 cents).
% For time points between the anchor points, the shift amount is calculated
% from linear intepolation. For time points after the last anchor point,
% the pitch shift value from the last anchor point is maintained.
%
% Therefore, the example below means "Exert no perturbation in the first
% second after the crossing of the intensity threshold; then linearly 
% ramp up the pitch-shift ratio from none to 100 cents in a period of 1
% sec; finally hold the 100-cent shift until the end of the supra-threshold
% event".
params.timeDomainPitchShiftSchedule = [0, 1.0; 1, 1.0; 1.150, 1.0595; 2.0, 1.0595; 2.15, 1.0];
params.rmsThresh = 0.011;

AudapterIO('init', params);
AudapterIO('reset');   % Reset;

if ~isOnline
    %% Offline demo: load and process input sound from wav file.
    
    mic = OfflineLoadBaselineVoice(dirs, expParam);

    micFrame = makecell(mic, expParam.frameLen);
    
    expParam.curSessTrial = [expParam.subject expParam.run];
    for n = 1 : length(micFrame)
        Audapter('runFrame', micFrame{n});
    end
else
    %% Online demo.
    Audapter('start');
    pause(4.5);
    Audapter('stop');
end

data = AudapterIO('getData');

dfSaveWavRec(data, expParam, dirs)

%% Visualize data.
if toDraw
    visPSSsigs(data)
end

% Play input and output sound.
if toPlay
    play_audio(data.signalIn, data.params.sr);
    play_audio(data.signalOut, data.params.sr);
end
end

function micproc = OfflineLoadBaselineVoice(dirs, expParam)

inputWav     = dirs.SavBaseFile;
[mic, micfs] = audioread(inputWav);

fs = expParam.sRate;

mic_resamp = resample(mic, fs, micfs);
micproc    = mic_resamp - mean(mic_resamp);
end

function visPSSsigs(data)

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0 : frameDur : frameDur * (size(data.fmts ,1) - 1);
% get the formant plot bounds

figure('Position', [200, 200, 800, 400]);
[s, f, t]    = spectrogram(data.signalIn, 64, 48, 1024, data.params.sr);
[s2, f2, t2] = spectrogram(data.signalOut, 64, 48, 1024, data.params.sr);

subplot(211);
imagesc(t, f, 10 * log10(abs(s))); 
axis xy;
colormap jet;
hold on;
plot(tAxis, data.fmts(:, 1 : 2), 'w', 'LineWidth', 2);
set(gca, 'YLim', [0, 4000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original');

subplot(212);
imagesc(t2, f2, 10 * log10(abs(s2)));
axis xy;
colormap jet;
hold on;
plot(tAxis, data.fmts(:, 1 : 2), 'w', 'LineWidth', 2);
plot(tAxis, data.sfmts(:, 1 : 2), 'g--', 'LineWidth', 2);
set(gca, 'YLim', [0, 4000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Shifted');

figure;
hold on;
sigTAxis1 = 0 : 1 / data.params.sr : ...
    1 / data.params.sr * (length(data.signalIn) - 1);
plot(sigTAxis1, data.signalIn, 'b-');
ylabel('Input waveform');
xlabel('Time (s)');

sigTAxis2 = 0 : 1 / data.params.sr : ...
    1 / data.params.sr * (length(data.signalOut) - 1);
plot(sigTAxis2, data.signalOut, 'r-');
legend({'Input', 'Output'});
ylabel('Output waveform');
xlabel('Time (s)');

figure;
hold on;
plot(sigTAxis1(1 : end - 1), diff(data.signalIn), 'b-');
ylabel('d Input waveform');
xlabel('Time (s)');
plot(sigTAxis2(1 : end - 1), diff(data.signalOut), 'r-');
legend({'Input', 'Output'});
ylabel('d Output waveform');
xlabel('Time (s)');

% When time-domain pitch shifting is activated, data.pitchHz and
% data.shfitedPitchHz record the tracked and shifted pitch values in Hz,
% for supra-threshold frames.
figure;
plot(tAxis, data.pitchHz, 'b-');
hold on;
plot(tAxis, data.shiftedPitchHz, 'r-');
legend({'Input', 'Output'});
xlabel('Time (s)');
ylabel('Pitch (Hz)');
end