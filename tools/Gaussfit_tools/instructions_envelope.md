### Working example for envelope measure ###

1. Go in the directory of the program envelope_measure.pro 

2. Type in a linux/mac terminal 'idl' to open an IDL console

3. Compile the program. Type (two times):
		 .compile envelope_measure
		 .compile envelope_measure

4. Check for errors during the compilation. If you had no errors, you should see the three following lines:
		% Compiled module: GET_NUMAX_AMAX.
		% Compiled module: ENVELOPE_MEASURE.
		% Compiled module: AUTOFILL.

5. Declare a string variable that contain the directory where the power spectrum of the stars are.
	For the moment, you have only Kepler-25 (KIC 4349452) in the 'Configurations/spec' directory.
	But you can put as many power spectrum as you have here (I am going to provide you these for the moment, at your request)
	For example, in my case, the directory is: '~/Dropbox/TOKYO/Programs/IDL/Configurations/spec/'
  Thus, in the IDL console (open in step 2) I type: 
  		dir_files='~/Dropbox/TOKYO/Programs/IDL/Configurations/spec/'
  		
6. Put the Identifier for the star(s) that you want to analyse, into the KIC_list variable. 
	In my case, I always use the KIC number for the identifier. But the identifier
	can be any kind of number (e.g. KOI number). Once you have decided for a convention, it is 
	better not to change it. 
   Thus in my case, I type the KIC number of Kepler-25:
   		KIC_list='4349452'
   		
7. Define the search range variable. For a main sequence you expect pulsation between 600 and 4000 microHz.
	You therefore have to type:
		search_range=[600, 4000] ; this is how we declare/initialize a 1d table in IDL

8. All input variables are now defined. The next step is to execute the main program, with 
	the proper arguments. Type;
		get_numax_Amax, dir_files, KIC_list, search_range

9. If all goes well, the last two lines in the IDL console will be something like:
		'Finished. A summary file -star_list.txt- was created in ~/Dropbox/TOKYO/Programs/IDL/Configurations/spec/'
		The program will exit now'

10. type 'exit' to leave the IDL console

11. Using the linux terminal, go to the directory where you spectrum is. In my example '~/Dropbox/TOKYO/Programs/IDL/Configurations/spec/'

12. by using the ls command, you should see two new files: 
		- star_list.txt
		- *4349452.sav.ps
	The star_list.txt has to be put into the '[...]/Preset-Analysis-v7.7/INPUT/' directory, to allow the next program to work...
    The .ps file is an image that could be open in you mac by using the 'open [name of the file to open]' command
    This .ps file allows you to check whether everything is alright regarding the way the program get_numax_Amax processed the data
    

Summary of the commands of the present example, assuming that you run a mac terminal:
	a. Mac terminal: cd ~/Dropbox/Programs/IDL/
	b. Linux terminal: idl
	c. IDL terminal: .compile envelope_measure
	d. IDL terminal: dir_files='~/Dropbox/Programs/IDL/Configurations/spec/'
	e. IDL terminal: KIC_list='4349452'
	f. IDL terminal: search_range=[600, 4000]
	g. IDL terminal: get_numax_Amax, dir_files, KIC_list, search_range
	h. IDL terminal: exit
	i. Mac terminal: cd ~/Dropbox/TOKYO/Programs/IDL/Configurations/spec/
	j. Mac terminal: ls
	k. Mac terminal: open 4349452.sav.ps
	
