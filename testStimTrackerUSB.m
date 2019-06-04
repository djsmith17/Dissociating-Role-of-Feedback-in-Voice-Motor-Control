clear all

s = serial('COM9');

fopen(s)

idn = fscanf(s, '%e', 14);

fclose(s)