###########################################################
# Part of ........
# Various functions
###########################################################

from Bio import PDB 
import numpy as np

######################################################################################################################

##########################
# Find_alignname_frompdb # 
##########################
# given PDB file and alignment, search the PDB for the UNIPROT name
# and then search the alignment for the complete name (with residues numbers) 
# calls Find_uniprot_frompdb AND Find_uniprot_completename
def Find_alignname_frompdb(align,pdb):
	uniname=Find_uniprot_frompdb(pdb) 
	if uniname == 'NOT_FOUND_IN_PDB':
		return uniname
	completename = Find_uniprot_completename(align,uniname);
	return completename
	

########################
# Find_uniprot_frompdb # 
########################
# Finds the uniprot name in the pdb file
def Find_uniprot_frompdb(pdb):
	with open(pdb,'r') as doc1:
		for line in doc1:
			line = line.split()
			if not line:
				continue
			if line[0]=='DBREF' and 'UNP' in line:
				uniname=line[line.index('UNP')+2]
				print '\tUNIPROT name found in PDB file: %s' % uniname
				return uniname
		print '\t!!ERROR!! UNIPROT name not found in PDB file'
		return 'NOT_FOUND_IN_PDB'



#############################
# Find_uniprot_completename # 
#############################
# given alignment and uniprot name, looks in the alignment for the complete uniprot name (with residues numbers)
def Find_uniprot_completename(align,uniname):
	flag=0 # notfound
	with open(align,'r') as doc2:
		for line in doc2:
			if line[0] == '>':
				if line.find(uniname)>0:
					unicomplete=line[1:-1]
					print '\tName found in the alignment: %s' % unicomplete
					flag+=1
	if flag > 1 :
		print '\t!!ERROR!! More than one sequence in the alignment with same UNIPROT name'
		return 'DUPLICATED_FOUND_IN_ALIGNMENT'
	if flag == 0 :
		print '\t!!ERROR!! No sequence in the alignment with that UNIPROT name'
		return 'NOT_FOUND_IN_ALIGNMENT'
	return unicomplete

######################################################################################################################

###################
# Create_alphabet #
###################
# creates an alphabet list, as specified from option: ('upper') or ('lower')
def Create_alphabet(option):
	if option=='upper' :
		return ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	elif option=='lower' :
		return ['a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y']
	else :
		return []

#########################
# Uppercase_transformer #
#########################
# it modifies the given string replacing all lowercase letters with uppercase
def Uppercase_transformer(s):
	lo=Create_alphabet('lower')
	up=Create_alphabet('upper')
	di=dict(zip(lo,up))
	for i in di.keys() :
		s = s.replace(i,di[i]) 
	return s

####################
# Letter_to_number #
####################
# given a letter and the lowercase treatment, returns the correspondent number
def Letter_to_number(letter,lowercase_treatment=0):
	number = -4
	if lowercase_treatment == 0:	## Tratto tutto allo stesso modo, come se fosse uppercase
		if letter in ['A','a']:		number = 0
		elif letter in ['C','c']:	number = 1
		elif letter in ['D','d']:	number = 2
		elif letter in ['E','e']:	number = 3
		elif letter in ['F','f']:	number = 4
		elif letter in ['G','g']:	number = 5
		elif letter in ['H','h']:	number = 6
		elif letter in ['I','i']:	number = 7
		elif letter in ['K','k']:	number = 8
		elif letter in ['L','l']:	number = 9
		elif letter in ['M','m']:	number = 10
		elif letter in ['N','n']:	number = 11
		elif letter in ['P','p']:	number = 12
		elif letter in ['Q','q']:	number = 13
		elif letter in ['R','r']:	number = 14
		elif letter in ['S','s']:	number = 15
		elif letter in ['T','t']:	number = 16
		elif letter in ['V','v']:	number = 17
		elif letter in ['W','w']:	number = 18
		elif letter in ['Y','y']:	number = 19
		elif letter in ['-','.']:	number = 20
		elif letter in ["B","Z","J","X","U","O","b","z","j","x","u","o"]:	## ambiguous residues
			number = -1


	elif lowercase_treatment == 1:	## Distinguo gli uppercase dai lowercase	
		if letter in ["a","r","n","d","c","e","q","g","h","i","l","k","m","f","p","s","t","w","y","v"]:	## lowercase residues
			number = -2
		elif letter == 'A':	number = 0
		elif letter == 'C':	number = 1
		elif letter == 'D':	number = 2
		elif letter == 'E':	number = 3
		elif letter == 'F':	number = 4
		elif letter == 'G':	number = 5
		elif letter == 'H':	number = 6
		elif letter == 'I':	number = 7
		elif letter == 'K':	number = 8
		elif letter == 'L':	number = 9
		elif letter == 'M':	number = 10
		elif letter == 'N':	number = 11
		elif letter == 'P':	number = 12
		elif letter == 'Q':	number = 13
		elif letter == 'R':	number = 14
		elif letter == 'S':	number = 15
		elif letter == 'T':	number = 16
		elif letter == 'V':	number = 17
		elif letter == 'W':	number = 18
		elif letter == 'Y':	number = 19
		elif letter == '-':	number = 20
		elif letter == ".":	## dot residue
			number = -3
		elif letter in ["B","Z","J","X","U","O","b","z","j","x","u","o"]:	## ambiguous residues
			number = -1
	
	if letter == "WT":	number = -10
	return number

####################
# Number_to_letter #
####################
# Given a number, it returns the letter or three-letter aa
def Number_to_letter(number,mod=1):
	if(mod==1):
		if number == 0:		letter = 'A'
		elif number == 1:	letter = 'C'
		elif number == 2:	letter = 'D'
		elif number == 3:	letter = 'E'
		elif number == 4:	letter = 'F'
		elif number == 5:	letter = 'G'
		elif number == 6:	letter = 'H'
		elif number == 7:	letter = 'I'
		elif number == 8:	letter = 'K'
		elif number == 9:	letter = 'L'
		elif number == 10:	letter = 'M'
		elif number == 11:	letter = 'N'
		elif number == 12:	letter = 'P'
		elif number == 13:	letter = 'Q'
		elif number == 14:	letter = 'R'
		elif number == 15:	letter = 'S'
		elif number == 16:	letter = 'T'
		elif number == 17:	letter = 'V'
		elif number == 18:	letter = 'W'
		elif number == 19:	letter = 'Y'	
		elif number == 20:	letter = '-'
		else :			letter = 'others'
	else:
		if number == 0:		letter = 'ALA'
		elif number == 1:	letter = 'CYS'
		elif number == 2:	letter = 'ASP'
		elif number == 3:	letter = 'GLU'
		elif number == 4:	letter = 'PHE'
		elif number == 5:	letter = 'GLY'
		elif number == 6:	letter = 'HIS'
		elif number == 7:	letter = 'ILE'
		elif number == 8:	letter = 'LYS'
		elif number == 9:	letter = 'LEU'
		elif number == 10:	letter = 'MET'
		elif number == 11:	letter = 'ASN'
		elif number == 12:	letter = 'PRO'
		elif number == 13:	letter = 'GLN'
		elif number == 14:	letter = 'ARG'
		elif number == 15:	letter = 'SER'
		elif number == 16:	letter = 'THR'
		elif number == 17:	letter = 'VAL'
		elif number == 18:	letter = 'TRP'
		elif number == 19:	letter = 'TYR'	
		elif number == 20:	letter = '---'
		else :			letter = 'others'
	return letter

################
# Three_to_one #
################
# given a 3-letters name of a residue, it returns the 1-letter identifier
def Three_to_one(residue3):
	if residue3 == 'ALA':	residue1 = 'A'
	elif residue3 == 'CYS':	residue1 = 'C'
	elif residue3 == 'ASP':	residue1 = 'D'
	elif residue3 == 'GLU':	residue1 = 'E'
	elif residue3 == 'PHE':	residue1 = 'F'
	elif residue3 == 'GLY':	residue1 = 'G'
	elif residue3 == 'HIS':	residue1 = 'H'
	elif residue3 == 'ILE':	residue1 = 'I'
	elif residue3 == 'LYS':	residue1 = 'K'
	elif residue3 == 'LEU':	residue1 = 'L'
	elif residue3 == 'MET':	residue1 = 'M'
	elif residue3 == 'ASN':	residue1 = 'N'
	elif residue3 == 'PRO':	residue1 = 'P'
	elif residue3 == 'GLN':	residue1 = 'Q'
	elif residue3 == 'ARG':	residue1 = 'R'
	elif residue3 == 'SER':	residue1 = 'S'
	elif residue3 == 'THR':	residue1 = 'T'
	elif residue3 == 'VAL':	residue1 = 'V'
	elif residue3 == 'TRP':	residue1 = 'W'
	elif residue3 == 'TYR':	residue1 = 'Y'
	else :
		residue1 = '-'
		#print residue3
		#sys.exit(-1)
	return residue1

################
# One_to_three #
################
# given a 1-letter name of a residue, it returns the 3-letters identifier
def One_to_three(residue1):
	if   residue1 == 'A': residue3 = 'ALA'
	elif residue1 == 'C': residue3 = 'CYS'	
	elif residue1 == 'D': residue3 = 'ASP'
	elif residue1 == 'E': residue3 = 'GLU'	
	elif residue1 == 'F': residue3 = 'PHE'	
	elif residue1 == 'G': residue3 = 'GLY'	
	elif residue1 == 'H': residue3 = 'HIS'	
	elif residue1 == 'I': residue3 = 'ILE'	
	elif residue1 == 'K': residue3 = 'LYS'	
	elif residue1 == 'L': residue3 = 'LEU'	
	elif residue1 == 'M': residue3 = 'MET'	
	elif residue1 == 'N': residue3 = 'ASN'	
	elif residue1 == 'P': residue3 = 'PRO'
	elif residue1 == 'Q': residue3 = 'GLN'
	elif residue1 == 'R': residue3 = 'ARG'	
	elif residue1 == 'S': residue3 = 'SER'	
	elif residue1 == 'T': residue3 = 'THR'	
	elif residue1 == 'V': residue3 = 'VAL'	
	elif residue1 == 'W': residue3 = 'TRP'	
	elif residue1 == 'Y': residue3 = 'TYR'
	else :
		residue3 = '---'
		#print residue3
		#sys.exit(-1)
	return residue3

#####################
# Matrix_to_numbers #
#####################
# Given a matrix with sequences in letters (MSA), returns a matrix in numbers
def Matrix_to_numbers(a,M,L,lowtre):
	na = np.zeros((M,L),dtype=int)
	for i in range(0,M):
	        for j in range(0,L):
	                na[i,j]=Letter_to_number(a[i,j],lowtre)
	return na	

################
# AA_to_colors #
################
# given a 3-letters name of a residue, it returns the 1-letter identifier
def AA_to_colors(residue3):
	if residue3 == 'ALA':	color = 'red'
	elif residue3 == 'CYS':	color = 'darkgreen'
	elif residue3 == 'ASP':	color = 'black'
	elif residue3 == 'GLU':	color = 'royalblue'
	elif residue3 == 'PHE':	color = 'aquamarine'
	elif residue3 == 'GLY':	color = 'yellow'
	elif residue3 == 'HIS':	color = 'khaki'
	elif residue3 == 'ILE':	color = 'fuchsia'
	elif residue3 == 'LYS':	color = 'darkviolet'
	elif residue3 == 'LEU':	color = 'silver'
	elif residue3 == 'MET':	color = 'navy'
	elif residue3 == 'ASN':	color = 'lime'
	elif residue3 == 'PRO':	color = 'white'
	elif residue3 == 'GLN':	color = 'dimgray'
	elif residue3 == 'ARG':	color = 'rosybrown'
	elif residue3 == 'SER':	color = 'pink'
	elif residue3 == 'THR':	color = 'maroon'
	elif residue3 == 'VAL':	color = 'navajowhite'
	elif residue3 == 'TRP':	color = 'steelblue'
	elif residue3 == 'TYR':	color = 'darkgoldenrod'
	return color

################
# Name_to_html #
################
# maps the name of the color to exadecimal number
def Name_to_html(nam):
	cnames = {
'aliceblue':            '#F0F8FF',
'antiquewhite':         '#FAEBD7',
'aqua':                 '#00FFFF',
'aquamarine':           '#7FFFD4',
'azure':                '#F0FFFF',
'beige':                '#F5F5DC',
'bisque':               '#FFE4C4',
'black':                '#000000',
'blanchedalmond':       '#FFEBCD',
'blue':                 '#0000FF',
'blueviolet':           '#8A2BE2',
'brown':                '#A52A2A',
'burlywood':            '#DEB887',
'cadetblue':            '#5F9EA0',
'chartreuse':           '#7FFF00',
'chocolate':            '#D2691E',
'coral':                '#FF7F50',
'cornflowerblue':       '#6495ED',
'cornsilk':             '#FFF8DC',
'crimson':              '#DC143C',
'cyan':                 '#00FFFF',
'darkblue':             '#00008B',
'darkcyan':             '#008B8B',
'darkgoldenrod':        '#B8860B',
'darkgray':             '#A9A9A9',
'darkgreen':            '#006400',
'darkkhaki':            '#BDB76B',
'darkmagenta':          '#8B008B',
'darkolivegreen':       '#556B2F',
'darkorange':           '#FF8C00',
'darkorchid':           '#9932CC',
'darkred':              '#8B0000',
'darksalmon':           '#E9967A',
'darkseagreen':         '#8FBC8F',
'darkslateblue':        '#483D8B',
'darkslategray':        '#2F4F4F',
'darkturquoise':        '#00CED1',
'darkviolet':           '#9400D3',
'deeppink':             '#FF1493',
'deepskyblue':          '#00BFFF',
'dimgray':              '#696969',
'dodgerblue':           '#1E90FF',
'firebrick':            '#B22222',
'floralwhite':          '#FFFAF0',
'forestgreen':          '#228B22',
'fuchsia':              '#FF00FF',
'gainsboro':            '#DCDCDC',
'ghostwhite':           '#F8F8FF',
'gold':                 '#FFD700',
'goldenrod':            '#DAA520',
'gray':                 '#808080',
'green':                '#008000',
'greenyellow':          '#ADFF2F',
'honeydew':             '#F0FFF0',
'hotpink':              '#FF69B4',
'indianred':            '#CD5C5C',
'indigo':               '#4B0082',
'ivory':                '#FFFFF0',
'khaki':                '#F0E68C',
'lavender':             '#E6E6FA',
'lavenderblush':        '#FFF0F5',
'lawngreen':            '#7CFC00',
'lemonchiffon':         '#FFFACD',
'lightblue':            '#ADD8E6',
'lightcoral':           '#F08080',
'lightcyan':            '#E0FFFF',
'lightgoldenrodyellow': '#FAFAD2',
'lightgreen':           '#90EE90',
'lightgray':            '#D3D3D3',
'lightpink':            '#FFB6C1',
'lightsalmon':          '#FFA07A',
'lightseagreen':        '#20B2AA',
'lightskyblue':         '#87CEFA',
'lightslategray':       '#778899',
'lightsteelblue':       '#B0C4DE',
'lightyellow':          '#FFFFE0',
'lime':                 '#00FF00',
'limegreen':            '#32CD32',
'linen':                '#FAF0E6',
'magenta':              '#FF00FF',
'maroon':               '#800000',
'mediumaquamarine':     '#66CDAA',
'mediumblue':           '#0000CD',
'mediumorchid':         '#BA55D3',
'mediumpurple':         '#9370DB',
'mediumseagreen':       '#3CB371',
'mediumslateblue':      '#7B68EE',
'mediumspringgreen':    '#00FA9A',
'mediumturquoise':      '#48D1CC',
'mediumvioletred':      '#C71585',
'midnightblue':         '#191970',
'mintcream':            '#F5FFFA',
'mistyrose':            '#FFE4E1',
'moccasin':             '#FFE4B5',
'navajowhite':          '#FFDEAD',
'navy':                 '#000080',
'oldlace':              '#FDF5E6',
'olive':                '#808000',
'olivedrab':            '#6B8E23',
'orange':               '#FFA500',
'orangered':            '#FF4500',
'orchid':               '#DA70D6',
'palegoldenrod':        '#EEE8AA',
'palegreen':            '#98FB98',
'paleturquoise':        '#AFEEEE',
'palevioletred':        '#DB7093',
'papayawhip':           '#FFEFD5',
'peachpuff':            '#FFDAB9',
'peru':                 '#CD853F',
'pink':                 '#FFC0CB',
'plum':                 '#DDA0DD',
'powderblue':           '#B0E0E6',
'purple':               '#800080',
'red':                  '#FF0000',
'rosybrown':            '#BC8F8F',
'royalblue':            '#4169E1',
'saddlebrown':          '#8B4513',
'salmon':               '#FA8072',
'sandybrown':           '#FAA460',
'seagreen':             '#2E8B57',
'seashell':             '#FFF5EE',
'sienna':               '#A0522D',
'silver':               '#C0C0C0',
'skyblue':              '#87CEEB',
'slateblue':            '#6A5ACD',
'slategray':            '#708090',
'snow':                 '#FFFAFA',
'springgreen':          '#00FF7F',
'steelblue':            '#4682B4',
'tan':                  '#D2B48C',
'teal':                 '#008080',
'thistle':              '#D8BFD8',
'tomato':               '#FF6347',
'turquoise':            '#40E0D0',
'violet':               '#EE82EE',
'wheat':                '#F5DEB3',
'white':                '#FFFFFF',
'whitesmoke':           '#F5F5F5',
'yellow':               '#FFFF00',
'yellowgreen':          '#9ACD32'
}
	return cnames[nam]
