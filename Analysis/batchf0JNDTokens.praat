############################
#  Resynthesizes all the sound files in the
#  specified directory to have flat pitch
#  of the specified frequency.  Files are
#  saved in a specified directory.
############################

form Resynthize files to have flat pitch
	text tokenDir
	sentence Sound_file_extension
        positive numTokens
	## positive Resynthesis_pitch_(Hz) 150
endform

#Here, you make a listing of all the sound files in a directory.
Create Strings as file list... list 'tokenDir$'/*'sound_file_extension$'

filename$ = Get string... 1

#A sound file is opened from the listing:
Read from file... 'tokenDir$'/'filename$'
sound_one$ = selected$ ("Sound")

start = Get start time
end = Get end time

To Manipulation... 0.01 75 600

resynthesis_pitch = 110
for iToken to numTokens

	# Create a new pitch tier with the flat pitch:

	select Sound 'sound_one$'

	Create PitchTier... 'sound_one$' start end
	Add point... start resynthesis_pitch
	Add point... end resynthesis_pitch

	# Combine and save the resulting file:

	select Manipulation 'sound_one$'
	plus PitchTier 'sound_one$'
	Replace pitch tier
	select Manipulation 'sound_one$'
	Get resynthesis (PSOLA)
	Write to WAV file... 'tokenDir$''iToken''filename$'

	##select Sound 'sound_one$'
	##plus Manipulation 'sound_one$'
	select PitchTier 'sound_one$'
	Remove

	##select Sound 'sound_one$'
	resynthesis_pitch = resynthesis_pitch + 20
endfor

select all
Remove
	