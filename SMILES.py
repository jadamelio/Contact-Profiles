#!/user/bin/python


#Jack D'Amelio January 14th 2017

#inputs:"
#TTD: 	#Creates all possible combiniations of scaffolds, linkers and building blocks, placing each SMILE string into a seperate text folder
	#Calculates/Indentifies "pure" structures, using only one type of building block
	#Use Open Babel to create molecule SVGs(Bool controlled)
	#[HOLD]For each pure structure, calculates similarity to each structure, creates set of 3 arrays, indexing: SMILE string, Open Babel SVGs(Bool controlled) and Similarity to pure moleucle
	#Creates some sort of graphic or table



import subprocess
import os
import shutil
from difflib import SequenceMatcher
from sys import argv
import stacking
from tkowalcz import molstruct
from scipy.spatial import distance
import math
import numpy
import subprocess
import Centroid
def generateSVG(smiPath):
	
	os.system("obabel -ismi " + smiPath + " -O " + smiPath[:-4] + ".svg ---errorlevel 1")#---errorlevel 5")

path = os.getcwd() + "/"

#shutil.rmtree(path +"Lib")
def generateXYZ(smiPath, ogMonomer):
	#print smiPath
	os.system("obabel -ismi " + smiPath + " -O " + smiPath[:-4] + ".xyz --gen3d ---errorlevel 1")
	grabAndSpin( smiPath[:-4] + ".xyz")
	
	#unBendAndCenter(smiPath[:-4] + ".xyz", ogMonomer)
	return smiPath[:-4] + ".xyz"

def spinXYZ(ogPath):
	#print ogPath
	#os.system("obabel -ismi " + ogPath + " -O " + ogPath[:-4] + ".xyz --gen3d ---errorlevel 1")

	ogFile = open(ogPath, 'r')
	name = []
	x = []
	y = []
	z = []
	waste = ogFile.next()
	waste = ogFile.next()

	while True:

		
		try: 
			readIn = ogFile.next().split()
			name.append(readIn[0])
			x.append(readIn[1])
			y.append(readIn[2])
			z.append(readIn[3])
		except StopIteration:
			break

	ogFile.close()

	newFile = open(ogPath[:-4] + ".in", 'w')
	newFile.write("$comment"+ '\n')
	newFile.write("$molecule"+ '\n')
	newFile.write("0 1"+ '\n')

	for i in range(len(name)):
		newFile.write(str(name[i]) + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i])  + '\n')
	newFile.write("$end"+ '\n')

	newFile.write("$rem"+ '\n')
	newFile.write("exchange		pbe"+ '\n')
	newFile.write("basis	3-21G"+ '\n')
	newFile.write("max_scf_cycles	0"+ '\n')
	newFile.write("$end"+ '\n')
	newFile.close()

	os.system("qchem " + ogPath[:-4] + ".in" + " "+ ogPath[:-4] + ".out >/dev/null 2>/dev/null")



	newMol = molstruct.Molecule()
	newMol.acquire_from_qchem_output(ogPath[:-4] + ".out")
	newFile = open(ogPath[:-4] + ".xyz", 'w')
	numAtom = newMol.num_atoms()
	newFile.write(str(numAtom) + "\n" + "\n")

	for a in range(newMol.num_atoms()):
		molspecs = "%-2s%15f%15f%15f\n" % (newMol.coordinates.atomic_symbols[a], \
	                   newMol.coordinates.cartesians[a][0], newMol.coordinates.cartesians[a][1], \
	                   newMol.coordinates.cartesians[a][2])
		newFile.write(molspecs)
		#print(molspecs)
	# newFile = ogPath[:-4] + ".xyz"
	# framework = stacking.reorienetMonomer(newFile)
	# name = framework[0]
	# x = framework[1]
	# y = framework[2]
	# z = framework[3]

	# layer1 = stacking.geomol(name,x,y,z)
	# array = []
	# array.append(layer1)
	newFile.close()

		

	# stacking.writeXYZ(array, ogPath[:-4] + ".xyz")

	return ogPath[:-4] + ".xyz"

def generateBlankXYZ(smiPath):
	#print "obabel -ismi " + smiPath + " -O " + smiPath[:-4] + ".xyz"+ " ---errorlevel 1"
	os.system("obabel -ismi " + smiPath + " -O " + smiPath[:-4] + ".xyz"+ " ---errorlevel 1")
	return smiPath[:-4] + ".xyz"

def generateSMI(xyzPath):
	os.system("obabel -ixyz " + xyzPath + " -O " + xyzPath[:-4] + ".smi ---errorlevel 1")
	return xyzPath[:-4] + ".smi"

def generateChemicals(scaffoldsFile, linkerFile, building_blockFile, outPutFile, smiLibPath):
	#print "GRAND OLD TIME"
	#print scaffoldsFile
	#print linkerFile
	#print linkerFile
	#print building_blockFile
	#print outPutFile
	#print smiLibPath
	os.system("java -jar SmiLib.jar -s "  + scaffoldsFile + " -l " + linkerFile + " -b " + building_blockFile + " -f " + outPutFile)
	#babel -ismi library.smi -osdf -x3 | ./babel -isdf -osmi
	'''
	lines = []
	file = open(outPutFile, 'r')
	while True:
		try:
			readIn = file.next().split()
			lines.append(readIn)
		except StopIteration:
			break
	file.close()
	file = open(outPutFile, 'w')
	for i in range(len(lines)):
		file.write(lines[i][1] + "\n")
	file.close()
	#print "EE"
	#print outPutFile
	os.system("babel -ismi" + outPutFile + " -osdf -x3| babel -isdf -O woo.smi ---errorlevel 1" )
	file = open(outPutFile, 'w')
	for i in range(len(lines)):
		file.write(lines[i][0] + " " + lines[i][1] + " \n")
	file.close()
	'''
	return outPutFile
def isolateChemicals(listPath, graphicToggle, xyzToggle, depositPath, ogMonomer):
	path = depositPath
	#os.mkdir(path +"Lib")
	 
	#Cuts each SMILE string into its own file
	structurePaths = []
	counter = 0
	#p#rint listPath
	if xyzToggle:
		print "		Gotta make those xyz files...Please wait..."
	readIn = open(listPath, 'r')
	while True:
		
		try:
			line = readIn.next()
			newFile = open(path + "/structure" + str(counter) + ".smi", 'w')
			newFile.write(line.split()[1])
			##print line.split()[1]
			structurePaths.append(path  + "/structure" + str(counter) + ".smi")
			newFile.close()
			
		except StopIteration:
			waste = 1
			break


		counter+=1

	if graphicToggle:
		for i in range(len(structurePaths)):

			##print(structurePaths[len(structurePaths)-1])
			generateSVG(structurePaths[i])
	if xyzToggle:
		for i in range(len(structurePaths)):

			##print(structurePaths[len(structurePaths)-1])
			try:
				generateXYZ(structurePaths[i], ogMonomer)
			except RuntimeError as re:
				print('Sorry but this structure was not able to be oriented '
          			'HERE: : {}'.format(re.args[0]))
				os.remove(structurePaths[i][:-4]+".svg")
				structurePaths.pop()
	
	return structurePaths
	print "		Done making those xyz files..."
'''
def main():

	listPath = argv[1]
	linkerFile = argv[2]
	scaffoldFile = argv[3]
	building_blockFile = argv[4]
	smiLibPath = argv[5]
	graphicToggle = False
	xyzToggle = False
	if len(argv) > 6:
		graphicToggle = bool(argv[6])
	if len(argv) > 7:
		xyzToggle = bool(argv[7])

	generateChemicals(scaffoldFile,linkerFile, building_blockFile, "output.txt",smiLibPath)
	isolateChemicals(listPath, graphicToggle, xyzToggle, os.getcwd())
	
'''

def unBendAndCenter(substitutedCurvedXYZ, straightUnsubsitutedOriginalMonomer):
	realCentroid = normalize(straightUnsubsitutedOriginalMonomer)
	badApple = qchemReiorient(substitutedCurvedXYZ)
	qchemCurve = grabAndSpin(badApple)
	normalize(qchemCurve, realCentroid)
def grabAndSpin(monomerPath, depth = 0):
	normalize(monomerPath)
	#grabAndSpinXY(monomerPath)
	if depth%6 == 0:
		escape = 0
		escape += grabAndSpinZY(monomerPath)
		escape += grabAndSpinXY(monomerPath)
		escape += grabAndSpinXZ(monomerPath)
	elif depth%6 == 1:
		escape = 0
		escape += grabAndSpinZY(monomerPath)		
		escape += grabAndSpinXZ(monomerPath)	
		escape += grabAndSpinXY(monomerPath)
	elif depth%6 == 2:
		escape = 0
		escape += grabAndSpinXY(monomerPath)
		escape += grabAndSpinZY(monomerPath)		
		escape += grabAndSpinXZ(monomerPath)

	elif depth%6 == 3:
		escape = 0
		escape += grabAndSpinXY(monomerPath)		
		escape += grabAndSpinXZ(monomerPath)
		escape += grabAndSpinZY(monomerPath)
	elif depth%6 == 4:	
		escape = 0
		escape += grabAndSpinXZ(monomerPath)
		escape += grabAndSpinZY(monomerPath)
		escape += grabAndSpinXY(monomerPath)
	

	elif depth%6 == 5:	
		escape = 0
		escape += grabAndSpinXZ(monomerPath)
		escape += grabAndSpinXY(monomerPath)
		escape += grabAndSpinZY(monomerPath)
		
	
	
	if escape == 3:
		return monomerPath
	else:
		grabAndSpin(monomerPath, depth+1)

def grabAndSpinZY(monomerPath):
	inital = open(monomerPath,'r')
	waste = inital.next()
	waste = inital.next()
	atmName = []
	atmX = []
	atmY = []
	atmZ = []
	#print(waste)
	while True:
		try:
			readIn = inital.next().split()
			#print(readIn)
			atmName.append(readIn[0])
			atmX.append(float(readIn[1]))
			atmY.append(float(readIn[2]))
			atmZ.append(float(readIn[3]))
		except StopIteration:
			#print("eeey")
			break

	


	mostDistantAtoms = [-1,-1]
	distanceToBeat = 0

	for i in range(len(atmName)):
		outerLoopAtom = (atmX[i], atmY[i], atmZ[i])
		for j in range(len(atmName)-i):
			innerLoopAtom = (atmX[i+j], atmY[i+j], atmZ[i+j])
			if distance.euclidean(outerLoopAtom, innerLoopAtom) > distanceToBeat:
				distanceToBeat = distance.euclidean(outerLoopAtom, innerLoopAtom)
				mostDistantAtoms = [i, i+j]
	#print atmName[mostDistantAtoms[0]], atmName[mostDistantAtoms[1]] 
	#print atmX[mostDistantAtoms[0]], atmX[mostDistantAtoms[1]] 
	#print atmY[mostDistantAtoms[0]], atmY[mostDistantAtoms[1]] 
	#print atmZ[mostDistantAtoms[0]], atmZ[mostDistantAtoms[1]] 
	if atmY[mostDistantAtoms[0]] != atmY[mostDistantAtoms[1]] :
		
		if atmZ[mostDistantAtoms[0]] > atmZ[mostDistantAtoms[1]] :
			#print "here:"
			adjacentLength = atmZ[mostDistantAtoms[0]] - atmZ[mostDistantAtoms[1]]
			oppositeLength = atmY[mostDistantAtoms[0]] - atmY[mostDistantAtoms[1]]
			angle =  -numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			#print adjacentLength
			#print oppositeLength
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				
				newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) + atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i])
				
				newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) - atmY[i]*numpy.sin(numpy.deg2rad(angle)))
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
			#print angle
		elif atmZ[mostDistantAtoms[0]] < atmZ[mostDistantAtoms[1]]:
			#print "there:"
			oppositeLength = atmY[mostDistantAtoms[1]] - atmY[mostDistantAtoms[0]]
			adjacentLength = atmZ[mostDistantAtoms[1]] - atmZ[mostDistantAtoms[0]]
			#print adjacentLength
			#print oppositeLength
			angle =  numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				
				newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) - atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i])
				newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) + atmY[i]*numpy.sin(numpy.deg2rad(angle)))
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
		#	print angle
		else:
			print "ALLS WELL???"
			return 1
	else:
		#print "ASDSey"
		angle =  90
		newName = []
		newX = []
		newY = []
		newZ = []
		for i in range(len(atmName)):
			newName.append(atmName[i])
			newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) - atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
			
			newX.append(atmX[i])
			
			newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) + atmY[i]*numpy.sin(numpy.deg2rad(angle)))
			
		newFile = open(monomerPath, 'w')
		newFile.write(str(len(newName))+ "\n")
		newFile.write("Comment"+ "\n")
		for i in range(len(atmName)):
			newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
		newFile.close()
	return 0
def grabAndSpinXY(monomerPath):
	inital = open(monomerPath,'r')
	waste = inital.next()
	waste = inital.next()
	atmName = []
	atmX = []
	atmY = []
	atmZ = []
	#print(waste)
	while True:
		try:
			readIn = inital.next().split()
			#print(readIn)
			atmName.append(readIn[0])
			atmX.append(float(readIn[1]))
			atmY.append(float(readIn[2]))
			atmZ.append(float(readIn[3]))
		except StopIteration:
			#print("eeey")
			break

	


	mostDistantAtoms = [-1,-1]
	distanceToBeat = 0

	for i in range(len(atmName)):
		outerLoopAtom = (atmX[i], atmY[i], atmZ[i])
		for j in range(len(atmName)-i):
			innerLoopAtom = (atmX[i+j], atmY[i+j], atmZ[i+j])
			if distance.euclidean(outerLoopAtom, innerLoopAtom) > distanceToBeat:
				distanceToBeat = distance.euclidean(outerLoopAtom, innerLoopAtom)
				mostDistantAtoms = [i, i+j]
	#print atmName[mostDistantAtoms[0]], atmName[mostDistantAtoms[1]] 
	#print atmX[mostDistantAtoms[0]], atmX[mostDistantAtoms[1]] 
	#print atmY[mostDistantAtoms[0]], atmY[mostDistantAtoms[1]] 
	#print atmZ[mostDistantAtoms[0]], atmZ[mostDistantAtoms[1]] 

	if atmY[mostDistantAtoms[0]] != atmY[mostDistantAtoms[1]] :
		
		if atmX[mostDistantAtoms[0]] > atmX[mostDistantAtoms[1]] :
			#print "here:"
			adjacentLength = atmY[mostDistantAtoms[0]] - atmY[mostDistantAtoms[1]]
			oppositeLength = atmX[mostDistantAtoms[0]] - atmX[mostDistantAtoms[1]]
			angle =  -numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			#print adjacentLength
			#print oppositeLength
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) - atmX[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) + atmY[i]*numpy.sin(numpy.deg2rad(angle)))
				newZ.append(atmZ[i])
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
			#print angle
		elif atmX[mostDistantAtoms[0]] < atmX[mostDistantAtoms[1]]:
			#print "there:"
			oppositeLength = atmX[mostDistantAtoms[1]] - atmX[mostDistantAtoms[0]]
			adjacentLength = atmY[mostDistantAtoms[1]] - atmY[mostDistantAtoms[0]]
			#print adjacentLength
			#print oppositeLength
			angle =  numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) + atmX[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) - atmY[i]*numpy.sin(numpy.deg2rad(angle)))
				newZ.append(atmZ[i])
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
		#	print angle
		else:
			print "ALLS WELL!!!"
			return 1
	else:
		#print "ASSDGGGy"
		angle =  90
		newName = []
		newX = []
		newY = []
		newZ = []
		for i in range(len(atmName)):
			newName.append(atmName[i])
			newY.append(atmY[i]*numpy.cos(numpy.deg2rad(angle)) + atmX[i]*numpy.sin(numpy.deg2rad(angle)))
			newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) - atmY[i]*numpy.sin(numpy.deg2rad(angle)))
			newZ.append(atmZ[i])
			newFile = open(monomerPath, 'w')
		newFile.write(str(len(newName))+ "\n")
		newFile.write("Comment"+ "\n")
		for i in range(len(atmName)):
			newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
		newFile.close()
	return 0
def grabAndSpinXZ(monomerPath):
	
	inital = open(monomerPath,'r')
	waste = inital.next()
	waste = inital.next()
	atmName = []
	atmX = []
	atmY = []
	atmZ = []
	#print(waste)
	while True:
		try:
			readIn = inital.next().split()
			#print(readIn)
			atmName.append(readIn[0])
			atmX.append(float(readIn[1]))
			atmY.append(float(readIn[2]))
			atmZ.append(float(readIn[3]))
		except StopIteration:
			#print("eeey")
			break

	


	mostDistantAtoms = [-1,-1]
	distanceToBeat = 0

	for i in range(len(atmName)):
		outerLoopAtom = (atmX[i], atmY[i], atmZ[i])
		for j in range(len(atmName)-i):
			innerLoopAtom = (atmX[i+j], atmY[i+j], atmZ[i+j])
			if distance.euclidean(outerLoopAtom, innerLoopAtom) > distanceToBeat:
				distanceToBeat = distance.euclidean(outerLoopAtom, innerLoopAtom)
				mostDistantAtoms = [i, i+j]


	centerY = (atmY[mostDistantAtoms[0]] + atmY[mostDistantAtoms[1]])/2
	cList = []
#	print centerY
#	print "______________"
	for i in range(len(atmName)):
		#print atmName[i] , atmY[i] 
		if atmName[i] == 'O': #and  atmY[i] - centerY < .5:
			cList.append(i)
	#print cList
	closestPair = []
	distanceToBeat = 100
	for i in range(len(cList)):
		outerLoopAtom = (atmX[cList[i]], atmY[cList[i]], atmZ[cList[i]])
		for j in range(len(cList)-i-1):

			innerLoopAtom = (atmX[cList[i+j+1]], atmY[cList[i+j+1]], atmZ[cList[i+j+1]])
			if distance.euclidean(outerLoopAtom, innerLoopAtom) < distanceToBeat:
				distanceToBeat = distance.euclidean(outerLoopAtom, innerLoopAtom)
				closestPair = [cList[i], cList[i+j+1]]

	#print atmName[closestPair[0]], atmName[closestPair[1]] 
	#print closestPair[0], closestPair[1]
	#print atmX[closestPair[0]], atmX[closestPair[1]] 
	#print atmY[closestPair[0]], atmY[closestPair[1]] 
	#print atmZ[closestPair[0]], atmZ[closestPair[1]] 
	#print "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
	#sprint closestPair
	if atmX[closestPair[0]] != atmX[closestPair[1]] :
		#print "EEE"
		if atmZ[closestPair[0]] > atmZ[closestPair[1]] :
			#print "here:"
			adjacentLength = atmZ[closestPair[0]] - atmZ[closestPair[1]]
			oppositeLength = atmX[closestPair[0]] - atmX[closestPair[1]]
			angle =  -numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			#print adjacentLength
			#print oppositeLength
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) - atmX[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) + atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
				newY.append(atmY[i])
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
			#print angle
		elif atmZ[closestPair[0]] < atmZ[closestPair[1]]:
			#print "there:"
			oppositeLength = atmX[closestPair[1]] - atmX[closestPair[0]]
			adjacentLength = atmZ[closestPair[1]] - atmZ[closestPair[0]]
			#print adjacentLength
			#print oppositeLength
			angle =  numpy.rad2deg(numpy.arctan(oppositeLength/adjacentLength))
			newName = []
			newX = []
			newY = []
			newZ = []
			for i in range(len(atmName)):
				newName.append(atmName[i])
				newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) - atmX[i]*numpy.sin(numpy.deg2rad(angle)))
				newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) + atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
				newY.append(atmY[i])
				newFile = open(monomerPath, 'w')
			newFile.write(str(len(newName))+ "\n")
			newFile.write("Comment"+ "\n")
			for i in range(len(atmName)):
				newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
			newFile.close()
		#	print angle
		else:
			#print "ALLS WELL###"
			return 1
	else:
		#print "print ey"
		angle =  90
		newName = []
		newX = []
		newY = []
		newZ = []
		for i in range(len(atmName)):
			newName.append(atmName[i])
			newZ.append(atmZ[i]*numpy.cos(numpy.deg2rad(angle)) + atmX[i]*numpy.sin(numpy.deg2rad(angle)))
			newX.append(atmX[i]*numpy.cos(numpy.deg2rad(angle)) - atmZ[i]*numpy.sin(numpy.deg2rad(angle)))
			newY.append(atmY[i])
			newFile = open(monomerPath, 'w')
		newFile.write(str(len(newName))+ "\n")
		newFile.write("Comment"+ "\n")
		for i in range(len(atmName)):
			newFile.write(str(newName[i]) + " "+  "{:10.4f}".format(newX[i])+ " " + "{:10.4f}".format(newY[i]) + " " + "{:10.4f}".format(newZ[i]) + "\n")
		newFile.close()
	return 0
def normalize(xyzPath, self = 'True'):
	oldFile = open(xyzPath, 'r')
	if self == 'True':
		centroid = stacking.centroid(xyzPath)
	else:
		centroid = self

	#print centroid
	name = []
	x = []
	y = []
	z = []
	waste = oldFile.next()
	waste = oldFile.next()
	while True:

		try:
			readIn = oldFile.next().split()
			name.append(readIn[0])

			x.append(float(readIn[1])-centroid[0])
			y.append(float(readIn[2])-centroid[1])
			z.append(float(readIn[3])-centroid[2])
			#print x
			#print y
			#print z
		except StopIteration:
			break

	oldFile.close
	oldFile = open(xyzPath, 'w')
	oldFile.write(str(len(name))+ " \n")
	oldFile.write("$comment"+ " \n")
	for i in range(len(name)):
		oldFile.write(str(name[i]) + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + " \n")
	oldFile.close()
	centroida = stacking.centroid(xyzPath)
#	print "final:"
	#print centroida
	return centroid

