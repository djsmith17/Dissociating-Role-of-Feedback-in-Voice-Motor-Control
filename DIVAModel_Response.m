
close all; clear all;

fs = 10; %samp/s

changeDown = 25:34;
lenChangeD = length(changeDown);
fpitch = 1/(2*lenChangeD);

t = 1:100;
time = ones(size(t));
pitch = ones(size(t));
pitch(changeDown) = 1 + -sin(2*pi*fpitch*(1:lenChangeD));

figure
plot(t, time, 'k')
hold on
plot(t, pitch, 'r')


axis([0 100 -5 5])