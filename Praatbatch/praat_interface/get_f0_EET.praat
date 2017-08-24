form get_f0
	text path 
	sentence Sound_file_extension 
	text outDir 
	positive voicing_threshold 
	positive pitch_floor 
	positive pitch_ceiling 
	positive max_formant 
	positive num_formant 
endform


Create Strings as file list... list 'path$'*'sound_file_extension$'
numberOfFiles = Get number of strings

# Go through all the sound files, one by one:

for ifile to numberOfFiles

	filename$ = Get string... ifile
	Read from file... 'path$''filename$'
	soundname$ = selected$ ("Sound", 1)
	tmin = Get start time
	tmax = Get end time
	To Pitch (ac)... 0.005 'pitch_floor' 15 off 0.03 'voicing_threshold' 0.01 0.35 0.14 'pitch_ceiling'
	select Sound 'soundname$'
	To Intensity... 100 0 yes
	select Sound 'soundname$'
	To Formant (burg)... 0.005 'num_formant' 'max_formant' 0.025 50

for i to (tmax-tmin)/0.005
    time = tmin + i * 0.005
    select Pitch 'soundname$'
    pitch = Get value at time... time Hertz Linear
    select Intensity 'soundname$'
    intensity = Get value at time... time Cubic
    select Formant 'soundname$'
    formant1 = Get value at time... 1 time Hertz Linear
    formant2 = Get value at time... 2 time Hertz Linear
    formant3 = Get value at time... 3 time Hertz Linear
    appendInfoLine: fixed$ (time, 2), " ", fixed$ (pitch, 3), " ", fixed$ (intensity, 3), " ", fixed$ (formant1, 3) , " ", fixed$ (formant2, 3) , " ", fixed$ (formant3, 3)
endfor
	
	
	titleline$ = "Time	F0(Hz)	Intensity(dB) F1(Hz) F2(Hz) F3(Hz) 'newline$'"
	fileappend 'outDir$'\'soundname$'.txt 'titleline$'
	fappendinfo 'outDir$'\'soundname$'.txt


	select Sound 'soundname$'
	plus Pitch 'soundname$'
        plus Formant 'soundname$'
        plus Intensity 'soundname$'
	Remove
	clearinfo
	select Strings list
endfor

#select all
#Remove

#fileappend "'outfile$'" 'titleline$'
	#fappendinfo 'outfile$' 