clear all

[SSN, fs] = audioread('SSN.wav');

SSN = repmat(SSN, 10, 1);

soundsc(SSN, fs)