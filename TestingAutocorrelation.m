%Program to find autocorrelation of a speech segment
clear; close all

[y, Fs] = audioread('Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Presentation\PrelimFiles\200Hz.wav'); %input: speech segment
y = y(1:480);

max_value = max(abs(y));
y = y/max_value;

t = (1/Fs:1/Fs:(length(y)/Fs))*1000;

sum1 = 0; 
autocor = 0;
for l = 0:(length(y)-1)
    sum1 = 0;
    for u = 1:(length(y)-l)
      s = y(u)*y(u+l);
      sum1 = sum1 + s;
    end
    autocor(l+1)= sum1;
end

kk = (1/Fs:1/Fs:(length(autocor)/Fs))*1000;

auto = autocor(21:160);
max1 = 0;
for uu = 1:140
    if(auto(uu) > max1)
      max1      = auto(uu);
      sample_no = uu;
    end 
end
pitch_period_To = (20+sample_no)*(1/Fs);
pitch_freq_Fo = 1/pitch_period_To;

figure
subplot(2,1,1);
plot(t,y);
xlabel('time in milliseconds')
title('A 30 millisecond segment of speech');

subplot(2,1,2);
plot(kk,autocor);
xlabel('time in milliseconds')
title('Autocorrelation of the 30 millisecond segment of speech');