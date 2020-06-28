file = uigetfile('*.wav', 'Select the wave file');

% Load the voice sample and noise
[data, fs] = audioread(file);
noise = randn(length(data),1);

% Normalize the voice sample
data  = data/sqrt((sum(data.^2))/length(data));

% Find the LPC of the voice sample
[a,g]    = lpc(data,6);
[Hv, Fv] = freqz(g, a, 256, fs);

% Filter and Normalize the noise
SSN   = filter(g, a, noise);
SSN   = SSN/sqrt((sum(SSN.^2))/length(data));

% Find the LPC of the noise
[a,g]    = lpc(SSN,6);
[Hn, Fn] = freqz(g, a, 256, fs);

figure
plot(Fv,log(abs(Hv)),'b');
hold on
plot(Fn,log(abs(Hn)),'r');