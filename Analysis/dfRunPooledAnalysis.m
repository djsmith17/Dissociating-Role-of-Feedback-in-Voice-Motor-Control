function dfRunPooledAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
pA.participants  = {'Pilot24'; 'Pilot24'}; %List of multiple participants.
pA.runs          = 'SF1';






end