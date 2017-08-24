# automate retrieval of pitch (f0)
# saves pulse timing data to f0.txt
# cooperates with MATLAB program get_f0.m
# this script can only run if PRAAT is running
# it cannot be run using PRAATCON
# 
# Joe Mendoza 2015
#

#setup the script
form get_f0
	sentence path
	sentence wav_fname
    	sentence outfile
	positive voicing_threshold
endform

clearinfo


Read from file... 'path$'

To Pitch (ac)... 0 75 15 off 0.03 'voicing_threshold' 0.01 0.35 0.14 600
Get mean... 0 0 Hertz

fappendinfo 'outfile$'

select all
Remove

clearinfo
