############################
#  Resynthesizes all the sound files in the
#  specified directory to have flat pitch
#  of the specified frequency.  Files are
#  saved in a specified directory.
############################

form Resynthize files to have flat pitch
	text tokenDir
	sentence Sound_file_extension
	text txtFileLoc
endform

#Here, you make a listing of all the sound files in a directory.
Create Strings as file list... list 'tokenDir$'/*'sound_file_extension$'

filename$ = Get string... 1

#A sound file is opened from the listing:
Read from file... 'tokenDir$'/'filename$'
sound_one$ = selected$ ("Sound")

start = Get start time
end = Get end time

To Pitch (ac)... 0.005 75 15 off 0.03 0.45 0.01 0.35 0.14 600

for i to (end - start)/0.005
    time = start + i * 0.005
    select Pitch 'sound_one$'
    pitch = Get value at time... time Hertz Linear
    appendInfoLine: fixed$ (time, 2), " ", fixed$ (pitch, 2)
endfor
	
	
#titleline$ = "Time	F0(Hz)	'newline$'"
#fileappend 'txtFileLoc$' 'titleline$'
fappendinfo 'txtFileLoc$'

select all
Remove

Quit