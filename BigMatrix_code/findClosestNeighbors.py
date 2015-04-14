import numpy as np
import csv
import timeit
import threading 


termToID = {}
IDtoTerm = {}
fileName = 	"Top100_closest_neighbors" + ".txt"
reportFile = open(fileName, "a") # append mode 

# to be filled out by functions of reading files
matrixU = np.matrix 
singularVals = []

'''
	============================================================
	Reading files code below 
'''


'''
	Make dictionaries mapping ID to terms, and vice versa 
	++ tested : YES
'''

def makeDictionaries():
	file_terms = "termList.txt"
	global termToID
	global IDtoTerm
	# build the dictionary for term list 		
	with open(file_terms, "rb") as f:
		reader = csv.reader(f, delimiter="\t")
		index = 0
		for row in reader:
			IDtoTerm[index] = row[0]
			termToID[row[0]] = index
			index += 1

	f.close()

'''
	Read in decomposed matrices U and V 
	++ tested : YES
'''
def readMatrix(filename, row, col): 
	global matrixU
	matrixU = np.zeros((row, col)) 
	with open(filename, "rb") as f:
		reader = csv.reader(f, delimiter = " ")
		x = []
		for row in reader:
			for i in row:
				if i != "":
					x.append(float(i))
			matrixU[x[0]] = x[1:] # x[0] is index and the rest are values 
			x = []
	f.close()


'''
	Read singular value files
	++ tested : YES
'''
def getSVals(filename, header = True):
	global singularVals
	with open(filename, "rb") as f:
		reader = csv.reader(f, delimiter = " ")
		index = 0
		for row in reader: 
			if (header == True and index == 0):
				index += 1 
			else:
				singularVals.append(float(row[0]))
	f.close()


''' 
	scale the matrix by square root such m_ij = sqrt(s_j)*m_ij
	for later's distance calculation to be sum(s_0*(a_0-b_0)^2 + s_1*(a_1-b_1)^2,...) 
''' 
def scaleMatrix():
	global matrixU
	matrixU = np.multiply(matrixU, np.sqrt(singularVals))


'''
 ===========================================================
 Analysis code below 
'''


'''
	@inputs: matrix - either U or V from SVD 
			 s - eigenvalues using to weight the score
			ID - the row index to compare to 
			 n - number of neighbors to find
			fun - distance function
			weighted - whether scores should be weighted by singular values
	@output: an array of n number of IDs 

'''
def findNeighbors(ID, n = 100):
	scores = []
	IDs = []
	n_rows = matrixU.shape[0]

	# get the score by comparing row by row and store as a tuple 
	#time1 = timeit.default_timer()
	for i in range(0, n_rows): 
		score = np.linalg.norm(matrixU[ID] - matrixU[i]) 
		scores.append((i, score))  # (drugID, score) tuple 
	#time2 = timeit.default_timer()

	# in place, sort by score smallest to largest 
	scores.sort(key = lambda tup: tup[1]) 
	#time3 = timeit.default_timer()

	#print "time to calculate all distances: " + str(time2-time1)
	#print "time to sort all the pairs: " + str(time3-time2)
	record(ID, scores, n)



'''
	@inputs: 
		ID of the term of interest 
		tuples: list of closest term and distance in the form of (ID, distance) 
		n: number of terms we want 
	@output: an array of n number of IDs 

'''
def record(ID, tuples, n):
	for i in range(0, n):
		reportFile.write(str(ID) + ' ' + str(tuples[i][0]) + ' ' + str(tuples[i][1]) + '\n')
	# in the form of (ID, ID of closest term, distance)

'''
	=======================================================================
'''



time1 = timeit.default_timer() # start time 
readMatrix("U_600", 163934, 529) 
time2 = timeit.default_timer()
print "To read in matrix :" + str(time2-time1)


makeDictionaries()
time3 = timeit.default_timer()
print "To build term dictionary: " + str(time3-time2)


getSVals("singularVals_600", True)
time4 = timeit.default_timer()
print "To read singular val: " + str(time4-time3)


scaleMatrix()
time5 = timeit.default_timer()
print "to scale the matrix: " + str(time5-time4)


for i in range(0, 163935):
	findNeighbors(i, 100)

time6 = timeit.default_timer()
print "To find top 100 terms: " + str(time6-time5)

print "Total run time: " + str(time5-time1)

reportFile.close()