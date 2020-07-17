import numpy as np

######################################
# Weight_calculator_forlongsequences # (by Cape)
######################################
# Freely inspired to calculate_evolutionary_constraints.m
# look above. It does not creates a Matrix, but calls each time a modified hamming function
# counts the weights, i.e. the number of similar sequences (similarity>threshold)
#cython: boundscheck=False
#cython: wraparound=False
def Weight_calculator_forlongsequences(alignment,threshold):
	M=alignment.shape[0]
	N=alignment.shape[1]
	cdef unsigned int i,j
	weight=np.zeros((M))
	cdef float inverse_threshold = (1. - threshold)*float(N)
	for i in range(0,M):
		for j in range(i+1,M):
			# Here I use the hamming function from scipy, without controls on vectors. (2x speedup)
			if (alignment[i,:] != alignment[j,:]).sum() <= inverse_threshold:
				weight[i] += 1
				weight[j] += 1
		weight[i] += 1
		print i,"over",M,"rows done...\r",
	return weight
