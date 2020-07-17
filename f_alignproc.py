###########################################################
# Part of ........
# Functions to modify the alignment
###########################################################

import numpy as np

############################
# Reduce_alignment_Columns #
############################
# returns a reduced alignment (matrix) starting from the number of the reference sequence
# it needs a matrix with NUMBERS not LETTERS
# reduces columns. 
def Reduce_alignment_Columns(ali,index):
	L_old=len(ali[0])
	results = Find_reduction_indexes_Columns(ali,index) # finds which columns to remove
	col_keep	= results[0]
	col_remove	= results[1]
	pfam_seq_num	= results[2]
	pfam_seq_index	= results[3]
	count_lower	= results[4]
	count_gap	= results[5]
	pfam_seq_num 	= map(None,pfam_seq_num)
	col_red_ali = ali[ :,col_keep ]	# reduction over columns	
	L_new = len(col_red_ali[0])
	print '\tColumns removed in the alignment: \t%d / %d (%.4f%%)' % ( len(col_remove),L_old,(len(col_remove)/float(L_old))*100 )
	print '\t\t%4d (%.4f%%) for gaps in the reference sequence' %( count_gap, (count_gap/float(L_old))*100)
	print '\t\t%4d (%.4f%%) for lowercase (if low treat is on)' %( count_lower, (count_lower/float(L_old))*100)
	return col_red_ali
#########################
# Reduce_alignment_Rows #
#########################
# returns a reduced alignment (matrix) starting from the number of the reference sequence
# it needs a matrix with NUMBERS not LETTERS
# reduces rows. 
def Reduce_alignment_Rows(ali,threshold):
	M_old=len(ali)
	results = Find_reduction_indexes_Rows(ali,threshold) # find which rows to remove
	row_keep	= results[0]
	row_remove	= results[1]
	count_lower	= results[2]
	count_ambig	= results[3]
	count_gap	= results[4]
	list_perc_gaps	= results[5]
	row_red_ali = ali[row_keep,:]	# reduction over rows	
	M_new = len(row_red_ali)
	print '\tRows (sequences) removed in the alignment: \t%d / %d (%.4f%%)' % ( len(row_remove),M_old,(len(row_remove)/float(M_old))*100 )
	print '\t\t%4d (%.4f%%) for gap percentage in the sequence above threshold' %(count_gap, (count_gap/float(M_old))*100)
	print '\t\t%4d (%.4f%%) for lowercase (if low treat is on)' %(count_lower, (count_lower/float(M_old))*100)
	print '\t\t%4d (%.4f%%) for ambiguous residues' %(count_ambig, (count_ambig/float(M_old))*100)
	return row_red_ali

##################################
# Find_reduction_indexes_Columns #
##################################
# finds which columns has to be removed from the aligment matrix
# out: [col_keep,col_remove,pfam_seq_num,pfam_seq_index,count_lower,count_gap]
def Find_reduction_indexes_Columns(ali,refer,treat=0):
	col_keep = np.zeros((0),dtype=int)	# indexes of columns to keep
	col_remove = np.zeros((0),dtype=int)	# indexes of columns to remove	
	pfam_seq_num = np.zeros((0),dtype=int)	
	pfam_seq_index = 0
	count_lower = 0
	count_gap = 0	
	for pos in xrange(0,len(ali[0])):
		if ali[refer,pos]==-2 :		# lowercase residue, remove  this has to be FIXED.
			col_remove = np.append(col_remove,pos)
			pfam_seq_num = np.append(pfam_seq_numb,-2)
			count_lower += 1
		elif ali[refer,pos] == -1 :	# ambiguous residue, remove
			print "RESIDUO AMBIGUO!!!"
			break
		elif ali[refer,pos] in [-3,20] :# gap residue '-' or '.', remove
			col_remove = np.append(col_remove,pos)
			count_gap +=1
		else:				# normal residue, keep
			col_keep = np.append(col_keep,pos)	
			if ali[refer,pos] >=0 :
				pfam_seq_num = np.append(pfam_seq_num,pfam_seq_index)
			pfam_seq_index += 1
	return [col_keep,col_remove,pfam_seq_num,pfam_seq_index,count_lower,count_gap]

###############################
# Find_reduction_indexes_Rows #
###############################
# finds which rows has to be removed from the aligment matrix
# out: [row_keep,row_remove,count_lower,count_ambig,count_gap,list_perc_gaps]
def Find_reduction_indexes_Rows(ali,threshold):
	L=len(ali[0])
	M=len(ali)	
	row_keep = np.zeros((0),dtype=int)		# indexes of rows to keep
	row_remove = np.zeros((0),dtype=int)		# indexes of rows to remove
	list_perc_gaps = np.zeros((0))
	count_lower = 0	
	count_ambig = 0
	count_gap = 0
	for sequence in xrange(0,M):
		verifica = 0
		for residue in ali[sequence,:]:
			if residue == -2:	# lowercase, break
				row_remove = np.append(row_remove,sequence)
				verifica = -1
				count_lower += 1
				break
			elif residue == -1:	# ambiguous, break
				row_remove = np.append(row_remove,sequence)
				verifica = -1
				count_ambig += 1
				break
		num_gaps = L - np.count_nonzero(ali[sequence,:]-20) + L - np.count_nonzero(ali[sequence,:]+3) # gaps are values 20 and -3
		perc =  num_gaps/float(L)	# percentage of gaps
		list_perc_gaps = np.append(list_perc_gaps,perc)	
		if (perc >= threshold and verifica==0) :		# if percentage is >= threshold, remove sequence
			row_remove = np.append(row_remove,sequence)
			count_gap += 1
			verifica = -1
		if verifica == 0 :
			row_keep = np.append(row_keep,sequence)
	return[row_keep,row_remove,count_lower,count_ambig,count_gap,list_perc_gaps]

