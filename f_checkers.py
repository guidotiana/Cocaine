###########################################################
# Part of ........
# Functions for the checking of parameters and errors
###########################################################
import sys
from f_printers import Print_usage
from f_printers import Print_pars

##############
# Check_args #
##############
# Control if all the arguments NEEDED to start the program have been specified
# if not, prints the usage returns FALSE. 
# if -help or -pars has been specified, it prints the necessary infos
def Check_args(a): 
	if a.argv[1]=='-pars': # if specified, it prints the list of all parameters
		Print_pars()
		return False
	if a.argv[1]=='-help': # if specified, it prints how to run the program
		Print_usage()
		return False
	if '-a' in a.argv and '-s' in a.argv and '-disorder' not in a.argv: # all right
		if '-p' not in a.argv:
			print '!WARNING! parameters file not specified (option -p), default values will be used.'
			print '\t for informations about parameters that can be specified from file, try ./main.py -pars'
		if '-uni' not in a.argv:
			print '!WARNING! UNIPROT name not specified (option -uni), automatic search will be tried (can lead to errors).'
		return True
	elif '-disorder' in a.argv and '-a' in a.argv and '-uni' in a.argv:
		print 'Disorder mode selected'
		if '-p' not in a.argv:
			print '!WARNING! parameters file not specified (option -p), default values will be used.'
			print '\t for informations about parameters that can be specified from file, try ./main.py -pars'
		if '-s' in a.argv:
			print '!WARNING! structure file specified (option -s). You have selected disorder mode (-disorder),'
			print 'thus the given structure will be ignored.'
		return True
	else:	# everything else, that is not correct
		Print_usage()
		return False			

#########################
# checkcreate_directory # 
#########################
# creates the directory for the output. If already present, it changes the name
def Checkcreate_directory(outdir):
	import os, shutil
	check = ""
	while os.path.exists(outdir) and check != 'y':
		print '!WARNING! Directory specified for output already exists:' 
		go = True
		while go:
			check = raw_input('do you want to overwrite it? please write y or n:\t')
			if check == 'y' or check == 'Y':
				shutil.rmtree(outdir)
				go = False
			if check == 'n' or check == 'N':
				outdir = raw_input('please specify a new output directory:\t')
				print 'Directory name changed to %s' % outdir
				go = False
	try: 
		os.makedirs(outdir)
		return outdir
	except OSError:
		if not os.path.isdir(outdir):
			raise NotImplementedError, 'Cannot find %s folder' % path


