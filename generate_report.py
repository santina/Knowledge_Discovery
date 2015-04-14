import numpy as np
import csv
from pprint import pformat

'''
	A line of the input file looks like:
	1163    1399    nadh    mental disorders        1
	The first and second columns are IDs for drug and disease (the text in the third and fourth columns).
	Then the fifth column is the count of sentences in which these two terms co-occur.
'''


# It's a sparse matrix with the below dimensions. The rows represent drugs, and the columns represent diseases.
rows=1563
cols=1521
filename='drug_disease_matrix.txt'

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
		IDtoDrug[int(row[0])] = row[2]    # dictionary
		IDtoDisease[int(row[1])] = row[3]  # dictionary
		DrugtoID[row[2]] = int(row[0])
		DiseasetoID[row[3]] = int(row[1])

fullMatrix = np.zeros((rows,cols))  # filling some zeros (placeholder content)
fullMatrix[xCoord,yCoord] = value   # making our matrix


# Then we do an SVD using numpy, false for reduced SVD 
rU, rs, rV = np.linalg.svd(fullMatrix, full_matrices=False) 
# print rU.shape, rs.shape, rV.shape  # (1563, 8) (8,) (8, 1521)

reconstructed = np.matrix

# try different distance functions. only Euclidean function can be weighted 
def calculateDistance(array1, array2, weights, fun = "Euclidean", weighted = False):
	sum = 0

	# measure Euclidean distance 
	if fun == "Euclidean": 
		for i in range(0, len(weights)):
			score = (array1[i]-array2[i])**2  # a contentious line 
			if weighted == True:
				score *= weights[i] # a contentious line 
			sum += score

	# measure cosine similarity -1: opposite, 1: exactly the same 
	elif fun == "cosine": 
		a1 = np.mat(array1)
		a2 = np.mat(array2)
		dotproduct = np.dot(a1, a2.T)
		cosine = dotproduct/(np.linalg.norm(a1)*np.linalg.norm(a2))
		sum += np.asarray(cosine)[0][0]

	return sum


''' 
	@input: drugsMatrix - the matrix that represents pattern of occurence for drug names 
	@input: drugID - the ID/index of the drug we wish to compare with all other drugs
	@input: weights - the eigenvalues or weights for the scores 
	@input: n - number of drugs we wish to return 
	@input: fun - distance function we're using 
	@output: an array of top n drugs ID for drugs with the most similar patterns of 
			occurence with drugID 
'''

def compareDrugs(drugsMatrix, drugID, weights, n, fun = "Euclidean", weighted = False):
	scores = []
	drugsIDs = []

	m = drugsMatrix.copy()
	n_rows = m.shape[0]

	# get the score by comparing row by row and store as a tuple 
	for i in range(0, n_rows): 
		score = calculateDistance(m[drugID], m[i], weights, fun, weighted) 
		scores.append((score, i))  # (score, drugID) tuple 
	
	# in place, sort by score smallest to largest 
	scores.sort(key = lambda tup: tup[0]) 

	# get the drugID of the top n 
	for i in range(0, n):
		drugsIDs.append(scores[i][1]) 

	return drugsIDs

''' 
	@input: diseaseMatrix - the matrix that represents pattern of occurence for disease names 
	@input: diseaseID - the ID/index of the disease we wish to compare with all other diseases
	@input: weights - an array of eigenvalues or weights for the scores 
	@input: n - number of disease we wish to return 
	@input: fun - distance function we're using 
	@output: an array of top n diseases ID for drugs with the most similar patterns of 
			occurence with diseaseID 
'''

def comparisonDisease(diseasesMatrix, diseaseID, weights, n, fun = "Euclidean", weighted = False):
	scores = []
	diseasesIDs = []
	m = diseasesMatrix.copy()
	m = np.matrix.transpose(m)
	
	n_rows = m.shape[0]

	# get the scores 
	for i in range(0, n_rows): 
		score = calculateDistance(m[diseaseID], m[i], weights, fun, weighted)
		scores.append((score, i))  # (score, drugID)
	
	scores.sort(key = lambda tup: tup[0]) # in place, in ascending order

	for i in range(0, n):
		diseasesIDs.append(scores[i][1]) # append the ID

	return diseasesIDs


''' 
	@input: drugArray - an array of drug IDs/indices
	@input: diseaseArray - an array of disease IDs/indices
	@output: a list of tuples representing each drug-disease pair that have positive values in 
		both the original co-occurence matrix and reconstructed matrix. Format of each item:  
		(drugName, diseaseName, score in the original matrix, score in the reconstructed matrix)
		sorted by score in the reconstructed matrix, in descending order 
'''
def findAssociation(drugArray, diseaseArray):
	tuples = []

	for i in drugArray:
		for j in diseaseArray:
			if fullMatrix[i][j] > 0 and reconstructed[i][j] > 0: 
				tuples.append((IDtoDrug[i], IDtoDisease[j], 
					fullMatrix[i][j], reconstructed[i][j]))

	tuples.sort(key = lambda tup: tup[3], reverse = True) 
	return tuples 

''' 
	@input: drugArray - an array of drug IDs/indices
	@input: diseaseArray - an array of disease IDs/indices
	@input: drugID - the ID/index of the drug of interest 
	@input: diseaseID - the ID/index of the disease of interest 
	@output: a list of tuples representing (specific drug , disease) and (drug, specific disease) pairs 
		that have positive values in both the original co-occurence matrix and reconstructed matrix. 
		Format of each item:  (drugName, diseaseName, score in the original matrix, score in the reconstructed matrix)
		sorted by score in the reconstructed matrix, in descending order 
'''

def findSpecificAssociation(drugArray, diseaseArray, drugID, diseaseID):
	tuples = []

	for i in drugArray:
		if fullMatrix[i][diseaseID] > 0 and reconstructed[i][diseaseID] > 0:
			tuples.append((IDtoDrug[i], IDtoDisease[diseaseID], 
				fullMatrix[i][diseaseID], reconstructed[i][diseaseID]))
			
	for j in diseaseArray:
		if fullMatrix[drugID][j] > 0 and reconstructed[drugID][j] > 0:
			tuples.append((IDtoDrug[drugID], IDtoDisease[j], 
				fullMatrix[drugID][j], reconstructed[drugID][j]))

	tuples.sort(key = lambda tup: tup[3], reverse = True) 
	return tuples 

''' 
	@input: drugID - the ID/index of the drug of interest 
	@input: diseaseID - the ID/index of the disease of interest 
	@input: distanceFunction -"Euclidean"(default), "cosine", 
	@input: n - number of drugs and disease to find 
	@input: toConsole - printing the report to console or not 
	@input: weighted - whether to weight the distance calculation
	@output: print statements and a report file 
'''
def generateReport(drugID, diseaseID, distFun = "Euclidean", n = 10, toConsole = False, weighted = True, numSV = 8):
	global rU, rV, rs, reconstructed

	# And then we extract only the top numSV vectors/values
	rU = rU[:,:numSV] # for drugs
	rV = rV[:numSV,:] # for diseases
	rs = rs[:numSV]  # this is already a array of eigenvalues 

	# Now we reconstruct the original matrix using this matrix decomposition we did above
	reconstructed = np.dot(np.dot(rU,np.diag(rs)), rV)

	similarDrugs = compareDrugs(rU, drugID, rs, n, distFun, weighted)
	similarDiseases = comparisonDisease(rV, diseaseID, rs, n, distFun, weighted)
	allpairs = findAssociation(similarDrugs, similarDiseases)
	specificPairs = findSpecificAssociation(similarDrugs, similarDiseases, drugID, diseaseID)
	
	header = "Inputs: \n" "drug: " + IDtoDrug[drugID] 
	header += "\ndisease: " + IDtoDisease[diseaseID] + "\n distance function: " + distFun
	header += "\n(score in the original matrix, score in reconstructed matrix) : \n" 
	header += "(" + str(fullMatrix[drugID][diseaseID]) + " , " + str(reconstructed[drugID][diseaseID]) + ")"
	header += "\nnumber of diseases and drugs to search: " + str(n) 
	header += "\nNumber of ranks: " + str(numSV)
	report = header 

	report += "\n======== Drugs that are similar to " + IDtoDrug[drugID] + "========\n"
	report += "drug IDs: \n" 
	report += pformat(similarDrugs) + "\n"
	report += "drug names: \n"
	for ID in similarDrugs: 
		report += IDtoDrug[ID] + ", "

	report += "\n\n======== Diseases that are similar to " + IDtoDisease[diseaseID] + "========\n"
	report += "Disease IDs: \n"
	report += pformat(similarDiseases) + "\n"
	report += "disease names: \n"
	for ID in similarDiseases:
		report += IDtoDisease[ID] + ", "

	report += "\n\n======== Associations of those similiar drugs and similar diseases =======\n"
	report += "(drug, disease, original score, reconstructed score)"
	report += "there are "+ str(len(allpairs)) + " pairs: \n"
	report += pformat(allpairs)
	

	report += "\n\n======== Associations containing our disease and drug of interest  ========\n"
	report += "there are "+ str(len(specificPairs)) + " pairs: \n"
	report += pformat(specificPairs)


	if (toConsole == True):
		print report 
	else:
		drugname = IDtoDrug[drugID].replace(" ", "")
		diseasename = IDtoDisease[diseaseID].replace(" ", "")
		filename = "reports/SVD_report_" + drugname + "_" + diseasename + str(numSV)+ ".txt"
		reportFile = open(filename, "w")
		reportFile.write(report)
		reportFile.close()


generateReport(DrugtoID["chloroquine"], DiseasetoID["anemia"], "Euclidean", 10, False, True, 8)
generateReport(DrugtoID["efalizumab"], DiseasetoID["spondylitis"], "Euclidean", 10, False, True, 42)
