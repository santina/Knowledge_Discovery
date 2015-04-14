

''' About this file 

	Given an input file that store the term occurance matrix w the format : drugID diseaseID drugName diseaseName 
	
	This python code chunk generate four files:
		closest_drugs.txt  : store top 100 closest drugs for each drug 
		closest_diseases.txt : store top 100 closest diseases for each disease 
		drug_IDs.txt : matching drug ID to drug terms 
		disease_IDs.txt : matching disease ID to disease term
'''


import numpy as np
import csv
import timeit

''' Loading the term matrix  

	A line of the input file looks like:
	1163    1399    nadh    mental disorders        1
		The first and second columns are IDs for drug and disease (the text in the third and fourth columns).
		Then the fifth column is the count of sentences in which these two terms co-occur.
'''

time0 = timeit.default_timer()

# The rows represent drugs, and the columns represent diseases.
rows=1563
cols=1521
filename='drug_disease_matrix.txt'

# This parameter controls the SVD later and has been set through testing on this data-set.
numSV = 8

# This code then loads all this data into a bunch of lists and dictionaries
xCoord = []
yCoord = []
value = []
IDtoDrug = {}
IDtoDisease = {}
DrugtoID = {}
DiseasetoID = {}

with open(filename, 'rb') as f:
	reader = csv.reader(f, delimiter="\t")
	for row in reader:
		# Extract the x and y coordinate for this count value
		xCoord.append(int(row[0])) # ID for drug
		yCoord.append(int(row[1])) # ID for disease
		value.append(int(row[4]))  # number of occurrences 
		
		# Store the drug name and disease too for later use, matching ID to drugs/diseases
		IDtoDrug[int(row[0])] = row[2]    
		IDtoDisease[int(row[1])] = row[3]  
		DrugtoID[row[2]] = int(row[0])
		DiseasetoID[row[3]] = int(row[1])
f.close()

fullMatrix = np.zeros((rows,cols))  # filling some zeros (placeholder content)
fullMatrix[xCoord,yCoord] = value   # making our matrix


# Then we do an SVD using numpy, false for reduced SVD 
rU, rs, rV = np.linalg.svd(fullMatrix, full_matrices=False) 


# And then we extract only the top numSV vectors/values
rU = rU[:,:numSV] # for drugs
rV = rV[:numSV,:] # for diseases
singularVals = rs[:numSV]   # an array of eigenvalues 
# print rU.shape, rs.shape, rV.shape  # (1563, 8) (8,) (8, 1521)

# Now we reconstruct the original matrix using this matrix decomposition we did above
reconstructed = np.dot(np.dot(rU,np.diag(singularVals)), rV)



''' 
	scale the matrix by square root such m_ij = sqrt(s_j)*m_ij
	for later's distance calculation to be sum(s_0*(a_0-b_0)^2 + s_1*(a_1-b_1)^2,...) 
''' 
def scaleMatrix():
	global matrixU
	matrixU = np.multiply(matrixU, np.sqrt(singularVals))

'''
	@inputs:
			ID - the row index to compare to 
			 n - number of neighbors to find
			fun - distance function
	@output: an array of n number of IDs 

'''
def findclosestDrugs(ID, n = 100, fun = "Euclidean"):
	scores = []
	IDs = []
	n_rows = rU.shape[0]

	# measure Euclidean distance 
	if fun == "Euclidean": 
		# scale the matrix first 
		matrixU = rU.copy()
		matrixU = np.multiply(matrixU, np.sqrt(singularVals))

		# Calculate distance from all other rows (drugs)
		for i in range(0, n_rows):
			score = np.linalg.norm(matrixU[ID] - matrixU[i]) 
			scores.append((i, score))  # (drugID, score) tuple 

	# measure cosine similarity -1: opposite, 1: exactly the same 
	elif fun == "cosine": 
		for i in range(0, n_rows):
			a1 = np.mat(rU[ID])
			a2 = np.mat(rU[i])
			dotproduct = np.dot(a1, a2.T)
			cosine = dotproduct/(np.linalg.norm(a1)*np.linalg.norm(a2))
			score = np.asarray(cosine)[0][0]
			scores.append((i, score))

	# in place, sort by score smallest to largest 
	scores.sort(key = lambda tup: tup[1]) 
	#time3 = timeit.default_timer()

	recordNeighbors(ID, scores, n, "MySQL_data/closest_drugs.txt")


def findclosestDiseases(ID, n = 100, fun = "Euclidean"):
	scores = []
	IDs = []
	n_cols = rV.shape[1]

	# measure Euclidean distance 
	if fun == "Euclidean": 
		# scale the matrix first 
		matrixV = rV.transpose().copy()
		matrixV = np.multiply(matrixV, np.sqrt(singularVals))

		# Calculate distance from all other rows (drugs)
		for i in range(0, n_cols):
			score = np.linalg.norm(matrixV[ID] - matrixV[i]) 
			scores.append((i, score))  # (drugID, score) tuple 

	# measure cosine similarity -1: opposite, 1: exactly the same 
	elif fun == "cosine": 
		for i in range(0, n_cols):
			a1 = np.mat(rV.transpose()[ID])
			a2 = np.mat(rV.transpose()[i])
			dotproduct = np.dot(a1, a2.T)
			cosine = dotproduct/(np.linalg.norm(a1)*np.linalg.norm(a2))
			score = np.asarray(cosine)[0][0]
			scores.append((i, score))

	# in place, sort by score smallest to largest 
	scores.sort(key = lambda tup: tup[1]) 
	#time3 = timeit.default_timer()

	recordNeighbors(ID, scores, n, "MySQL_data/closest_diseases.txt")



'''
	@inputs: 
		ID of the term of interest 
		tuples: list of closest term and distance in the form of (ID, distance) 
		n: number of terms we want 
		filename : the filename to write the result to 
	@output: an array of n number of IDs 

'''
def recordNeighbors(ID, tuples, n, filename):
	f = open(filename, "a")
	for i in range(0, n):
		f.write(str(ID) + '\t' + str(tuples[i][0]) + '\t' + str(tuples[i][1]) + '\t' + str(i) + '\n')
	# in the form of (ID, ID of closest term, distance, rank)


'''
	Create the ID to noun phrase file
'''
def writeDictionary():
	# write to disease_IDs.txt
	disease_IDs = open("MySQL_data/disease_IDs.txt", "a")
	for i in range(0, cols):
		disease_IDs.write(str(i) + '\t' + IDtoDisease[i] + '\n')
	disease_IDs.close()

	# write to drug_IDs.txt
	drug_IDs = open("MySQL_data/drug_IDs.txt", "a")
	for i in range(0, rows):
		drug_IDs.write(str(i) + '\t' + IDtoDrug[i] + '\n')
	drug_IDs.close()

time1 = timeit.default_timer()

# create dictionary 
writeDictionary()
# create top 100 closest neighbor terms 

for i in range(0, cols):
	findclosestDiseases(i)
for i in range(0, rows):
	findclosestDrugs(i)

time2 = timeit.default_timer()

print str(time1 - time0) + ' ' +  str(time2 - time1) + ' ' + str(time2 - time0)