# automate retrieval of pitch (f0)
# saves pulse timing data to f0.txt
# cooperates with MATLAB program get_f0.m
# this script can only run if PRAAT is running
# it cannot be run using PRAATCON
# 
# Joe Mendoza 2015
# Edited by Ashling Lupiani 2016

#setup the script
form get_f0
	sentence path
	sentence wav_fname
    	sentence outfile
	positive voicing_threshold
	positive pitch_floor
	positive pitch_ceiling
endform

clearinfo


Read from file... 'path$'

#For an explanation of the variables on the following line, see 'pitch parameter standards for praat.txt'
To Pitch (ac)... 0 'pitch_floor' 15 off 0.03 'voicing_threshold' 0.01 0.35 0.14 'pitch_ceiling'

Get mean... 0 0 Hertz

fappendinfo 'outfile$'

select all
Remove

clearinfo
