############################
#  Resynthesizes all the sound files in the
#  specified directory to have flat pitch
#  of the specified frequency.  Files are
#  saved in a specified directory.
############################

form Resynthize files to have flat pitch
	text baseTokenFile
	sentence Sound_file_extension
	text txtFileLoc
endform

#A sound file is opened from the listing:
Read from file... 'baseTokenFile$'
sound_one$ = selected$ ("Sound")

start = Get start time
end   = Get end time

To Pitch (ac)... 0.005 75 15 off 0.03 0.45 0.01 0.35 0.14 600

for i to (end - start)/0.005
    time = start + i * 0.005
    select Pitch 'sound_one$'
    pitch = Get value at time... time Hertz Linear
    appendInfoLine: fixed$ (time, 2), " ", fixed$ (pitch, 2)
endfor
	
fappendinfo 'txtFileLoc$'

select all
Remove

Quit