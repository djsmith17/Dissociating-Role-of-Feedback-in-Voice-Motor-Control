############################
#
#  Resynthesizes all the sound files in the
#  specified directory to have flat pitch
#  of the specified frequency.  Files are
#  saved in a specified directory.
#
#
# This is a version copied straight from the web. I am now going to make a new one and modify
############################


form Resynthize files to have flat pitch
	comment Directory of sound files
	text sound_directory C:\users\djsmith\Desktop
	sentence Sound_file_extension .wav
	comment Directory of finished files
	text end_directory C:\users\djsmith\Desktop\PraatSynthFiles\
	comment Resynthesize to this frequency
	positive Resynthesis_pitch_(Hz) 150
endform

# Here, you make a listing of all the sound files in a directory.

Create Strings as file list... list 'sound_directory$'/*'sound_file_extension$'
numberOfFiles = Get number of strings


for ifile to numberOfFiles
	filename$ = Get string... ifile

	# A sound file is opened from the listing:

	Read from file... 'sound_directory$'/'filename$'
	sound_one$ = selected$ ("Sound")

	To Manipulation... 0.01 60 400

	# Create a new pitch tier with the flat pitch:

	select Sound 'sound_one$'
	start = Get start time
	end = Get end time
	Create PitchTier... 'sound_one$' start end
	Add point... start resynthesis_pitch
	Add point... end resynthesis_pitch

	# Combine and save the resulting file:

	select Manipulation 'sound_one$'
	plus PitchTier 'sound_one$'
	Replace pitch tier
	select Manipulation 'sound_one$'
	Get resynthesis (PSOLA)
	Write to WAV file... 'end_directory$''filename$'

	select Sound 'sound_one$'
	plus Manipulation 'sound_one$'
	plus PitchTier 'sound_one$'
	Remove

	select Strings list
endfor

select all
Remove
	