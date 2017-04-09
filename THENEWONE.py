#!/user/bin/python

from damelij import SMILES
from damelij import stacking
from sys import argv
import numpy
import os
import matplotlib
import matplotlib.pyplot
import shutil
import pylab







def generateScaffoldFiles(bulkScaffoldPath):
	file = open(bulkScaffoldPath, 'r')
	scaffoldStrings = []
	while True:
		try:
			readIn = file.next().split()
			##print readIn
			scaffoldStrings.append(readIn[0])
		except StopIteration:
			break

	folderPath = os.getcwd() + "/scaffoldStorage"
	try:
		os.mkdir(folderPath)

	except OSError:
		pass
	newFilePaths = []
	for i in range(len(scaffoldStrings)):
		newFile = open(folderPath + "/scaffold" + str(i) + ".smi", 'w')
		newFile.write(scaffoldStrings[i])
		newFilePaths.append(folderPath + "/scaffold" + str(i) + ".smi")

	return newFilePaths

	
def generateBaseline(originalPath,scaffoldPath, vdwDataPath, blastCap, xOffset, yOffset, numSheests, savePath, mainCentroid):
	'''
	temp = "tempBaseline/"
	os.mkdir(temp)
	readIn = open(scaffoldPath, 'r')
	line = readIn.next()
	limit = len(line)
	killZone = []
	print "		Killing [R-] Sites..."
	for i in range(limit):
		if line[i] == '[' and line[i+1] == 'R'and line[i+3] == ']':
			killZone.append(i)
			killZone.append(i+1)
			killZone.append(i+2) 
			killZone.append(i+3)
		elif line[i] == 'F':
			killZone.append(i)	
	#print killZone
	newLine = line


	for i in range(len(killZone)):
		#print newLine
		newLine = newLine[:killZone[i]] + newLine[killZone[i]+1:]
		for j in range(len(killZone)):
			killZone[j] -= 1
			
	killZone = []
	for i in range(len(newLine)):
		if newLine[i] == '(' and newLine[i+1] == ')':
			killZone.append(i)
			killZone.append(i+1)
	evenNewerLine = newLine
	for i in range(len(killZone)):
		#print evenNewerLine	
		evenNewerLine = evenNewerLine[:killZone[i]] + evenNewerLine[killZone[i]+1:]
		for j in range(len(killZone)):
			killZone[j] -= 1
	#print evenNewerLine		


	newFile = open('tempBaseline/temp.smi', 'w')
	newFile.write(evenNewerLine)
	newFile.close()

	print "		Generating XYZ..."
	xyzPath = SMILES.generateXYZ('tempBaseline/temp.smi', mainCentroid)
	#xyzPath = SMILES.theRMSDTrick(xyzPath, originalPath)
	##print xyzPath + "
	'''#            OOOO"

	
	
	
	print "		Graphing..."
	#baseGraph = stacking.graphAnglebyOffset(originalPath, vdwDataPath, blastCap, xOffset, yOffset, numSheests)
	h =stacking.precisionGraph(originalPath, vdwDataPath, .1, 0,0,2)
	##print "OOOOOOOOOOOOOOO    " + savePath
	print "		Histographing... :D"
	#graphScatter(baseGraph, savePath + "baseProfile.png")
	l =coordsToText(h, savePath +"P.baseProfile.txt")
	#arrayGraphToHist(baseGraph, savePath + "baseProfile.png")
	print "		Textographing... :D"
	#graphToText(baseGraph, savePath +" baseProfile.txt")
	analyticalGraphFromFileToScatter(l, savePath + "P.baseProfile.png", "blue")
	#shutil.rmtree(temp)
	#return baseGraph


def getMinMaxList(xList,yList):

	#removes all plateaus, replaces with single point
	blastIndex = []
	counter  = 0
	while True:
		#print len(yList)
		#print len(xList)
		label = 0
		for i in range(1, len(xList)-1):
			#print i, len(xList)
			if i + 1 == len(xList)-1:
				break
			try: 
				#print i
				if yList[i] == yList[i+1] and yList[i] == yList[i-1]:
					#print '---'
					#print i,  len(yList)
					#print len(xList)
					#print "hi"
					xList.pop(i-1)
					#print "looo", xList[i+1]
					xList.pop(i+1)
					yList.pop(i-1)
					yList.pop(i+1)
					#print "////"
					counter += 1
					label = 1

					#print len(yList)
					#print len(xList)
					#print '---'
					#i == i - 2
				elif yList[i] == yList[i+1]:
					#print "ho"
					yList[i] = (yList[i]+yList[i+1])/2
					xList[i] = (xList[i]+xList[i+1])/2
					xList.pop(i+1)
					yList.pop(i+1)
					counter += 1
					label = 1
				#	print len(yList)
					#print len(xList)
					#i  == i - 1
				elif yList[i] == yList[i-1]:
					#print "lets go"
					yList[i] = (yList[i]+yList[i-1])/2
					xList[i] = (xList[i]+xList[i-1])/2
					xList.pop(i-1)
					yList.pop(i-1)
					counter += 1
					label = 1
				#	print len(yList)
				#	print len(xList)
					#i  == i - 1
			except IndexError:
				#break
				pass
		
		if label == 0:
			break
	#print len(yList)
	#print len(xList)
	while True:
		if yList[0] == yList[len(yList)-1]:
				yList[0] =(yList[0]+yList[len(yList)-1])/2
				xList[0] = -(-xList[0]+xList[len(yList)-1])/2
				yList.pop(len(yList)-1)
				xList.pop(len(yList)-1)	
				counter +=1	
		else:
			break
	print "Removed " + str(counter) + "data points!"
	#print len(xList)
	#print len(yList)
	minMaxListA = []
	minMaxListD = []
	identity = []
	#print xList

	#print yList[len(yList)-1], yList[0],  yList[1] 
	if yList[0] < yList[1] and yList[0] < yList[len(yList)-1]:
		#print "hi"
		minMaxListD.append(yList[0])

		minMaxListA.append(xList[0])
		identity.append(1)
	elif yList[0] > yList[1] and yList[0] > yList[len(yList)-1]:
		minMaxListD.append(yList[0])
		minMaxListA.append(xList[0])
		identity.append(-1)

	for i in range(1, len(xList)-1):
		if yList[i] < yList[i+1] and yList[i] < yList[i-1]:
			minMaxListD.append(yList[i])
			minMaxListA.append(xList[i])
			identity.append(1)
		elif yList[i] > yList[i+1] and yList[i] > yList[i-1]:
			minMaxListD.append(yList[i])
			minMaxListA.append(xList[i])
			identity.append(-1)


	if yList[len(yList)-1] < yList[0] and yList[len(yList)-1] < yList[len(yList)-2]:
		minMaxListD.append(yList[len(yList)-1])
		minMaxListA.append(xList[len(yList)-1])
		identity.append(1)
	elif yList[len(yList)-1] > yList[0] and yList[len(yList)-1] > yList[len(yList)-2]:
		minMaxListD.append(yList[len(yList)-1])
		minMaxListA.append(xList[len(yList)-1])
		identity.append(-1)
	#print yList
	#print yList[len(yList)-2] , yList[len(yList)-1],yList[0],yList[1]
	#print identity
	matplotlib.pyplot.scatter(minMaxListA,minMaxListD, color = "blue")
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-')
	matplotlib.pyplot.yticks(numpy.arange(3, 7+1, .5))
	matplotlib.pyplot.xticks(numpy.arange(-180, 180+1, 30))
	#matplotlib.pyplot.show()
	return [minMaxListA, minMaxListD, identity]

def subtractGraph(baseLineGraph, subtractingGraph):
	
	for i in range(len(subtractingGraph)):
	
		for j in range(len(subtractingGraph[i])):
			#print subtractingGraph[i][j], baseLineGraph[i][j], subtractingGraph[i][j]-baseLineGraph[i][j]
			print subtractingGraph[i][j]
			print baseLineGraph[i][j]
			print 'OOO'
			subtractingGraph[i][j]-=baseLineGraph[i][j]
			print subtractingGraph[i][j]
	return subtractingGraph

def normalizeGraph(newGraph):
	minV = 0
	for i in range(len(newGraph)):
		for j in range(len(newGraph[i])):
			if newGraph[i][j] < minV:
				minV = newGraph[i][j]
	for i in range(len(newGraph)):
		for j in range(len(newGraph[i])):
			newGraph[i][j] += -1*minV
	return newGraph
def arrayGraphToHist(graph,savePath):
	arrayGraphToHistColor(graph,"color-"+savePath)
	arrayGraphToHistColorBlur(graph,"colorblur-"+savePath)
	arrayGraphToHistGrey(graph,"grey-"+savePath)
	
def graphToText(graph, savePath):
	
	newFile = open(savePath, 'w')
	for i in range(len(graph)):
		for j in range(len(graph[i])):
			if graph[i][j] > 0:
				newFile.write(str(i) + " " + str(4 + (j*.25)) + " "  + str(graph[i][j]) + "\n")
	newFile.close()
	return savePath
def coordsToText(coords, savePath):
	
	newFile = open(savePath, 'w')
	for i in range(len(coords)):
		newFile.write(str(coords[i][0]) + " "  + str(coords[i][1]) + " 1  "+"\n")
	newFile.close()
	return savePath


def arrayGraphToHistColorBlur(graph, savePath):
	a = numpy.matrix(graph)
	matplotlib.pyplot.imshow(a, aspect = 'auto')
	 #,extent=[0,180,4,10], aspect = 'auto')
	#(a, interpolation = 'none', extent=[0,180,4,10])#, aspect = 'auto')
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	#matplotlib.pyplot.xticks(numpy.arange(4, 10+1, .25))
	#matplotlib.pyplot.yticks(numpy.arange(0, 360+1, 10))

	matplotlib.pyplot.savefig(savePath)
	matplotlib.pyplot.clear()
def arrayGraphToHistColor(graph, savePath):
	a = numpy.matrix(graph)
	matplotlib.pyplot.imshow(a, aspect = 'auto',interpolation='nearest')
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	#matplotlib.pyplot.xticks(numpy.arange(4, 10+1, .25))
	#matplotlib.pyplot.yticks(numpy.arange(0, 360+1, 10))
	 #,extent=[0,180,4,10], aspect = 'auto')
	#(a, interpolation = 'none', extent=[0,180,4,10])#, aspect = 'auto')
	pylab.xlim([4,10])
	
	matplotlib.pyplot.savefig(savePath)
	matplotlib.pyplot.clear()
def arrayGraphToHistGrey(graphA, savePath):
	graph = [row[:] for row in graphA]

	for i in range(len(graph)):
		for j in range(len(graph[i])):
			if graph[i][j] >= 1:
				graph[i][j] = 1
			else:
				graph[i][j] = 0
	a = numpy.matrix(graph)
	matplotlib.pyplot.imshow(a, cmap='Greys' ,aspect = 'auto',interpolation='nearest')
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	#matplotlib.pyplot.xticks(numpy.arange(4, 10+1, .25))
	#matplotlib.pyplot.yticks(numpy.arange(0, 360+1, 10))
	 #,extent=[0,180,4,10], aspect = 'auto')
	#(a, interpolation = 'none', extent=[0,180,4,10])#, aspect = 'auto')
	pylab.xlim([4,10])
	#pylab.lim([0,360])
	matplotlib.pyplot.savefig(savePath)
	matplotlib.pyplot.clear()

	#matplotlib.pyplot.matshow(a)
	#matplotlib.pyplot.show()

def graphFromFile(path):
	graph = [[0 for i in range(24)] for j in range(360)]
	file = open(path, 'r')
	while True:
		try:
			readIn = file.next().split()
			#print readIn

			graph[int(readIn[0])][int((float(readIn[1])-4)/.25)] = float(readIn[2])
		except StopIteration:
			break
	file.close()

	return graph

def reverseGraphFromFileToScatter(path, savePath):
	angle = []
	dist = []
	file = open(path, 'r')
	while True:
		try:
			readIn = file.next().split()
			#print readIn

			if readIn[2] > 0:
				angle.append(360-float(readIn[0]))
				dist.append(float(readIn[1]))
		except StopIteration:
			break
	file.close()
	matplotlib.pyplot.scatter(angle,dist)
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	matplotlib.pyplot.yticks(numpy.arange(4, 7+1, .1))
	matplotlib.pyplot.xticks(numpy.arange(0, 360+1, 5))
	#matplotlib.pyplot.show()
	matplotlib.pyplot.savefig(savePath)
	#matplotlib.pyplot.
	matplotlib.pyplot.clf()
	matplotlib.pyplot.cla()
	matplotlib.pyplot.close()
	#matplotlib.pyplot.show()

def graphFromFileToScatter(path, savePath):
	angle = []
	dist = []
	file = open(path, 'r')
	while True:
		try:
			readIn = file.next().split()
			#print readIn

			if readIn[2] > 0:
				angle.append(float(readIn[0]))
				dist.append(float(readIn[1]))
		except StopIteration:
			break
	file.close()
	matplotlib.pyplot.scatter(angle,dist)
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	matplotlib.pyplot.yticks(numpy.arange(0, 7+1, .5))
	matplotlib.pyplot.xticks(numpy.arange(-180, 180+1, 20))
	#matplotlib.pyplot.show()
	matplotlib.pyplot.savefig(savePath)
	#matplotlib.pyplot.
	matplotlib.pyplot.clf()
	matplotlib.pyplot.cla()
	matplotlib.pyplot.close()
	#matplotlib.pyplot.show()
	
def graphScatter(graph, savePath):
	angle = []
	dist = []

	for i in range(len(graph)):
		for j in range(len(graph[i])):
			if graph[i][j] > 0:
				angle.append(i)
				dist.append(4+(.25*j))
	matplotlib.pyplot.scatter(angle,dist)
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-', color = 'b')
	matplotlib.pyplot.yticks(numpy.arange(0, 10+1, .5))
	matplotlib.pyplot.xticks(numpy.arange(0, 360+1, 50))
	#matplotlib.pyplot.show()
	matplotlib.pyplot.savefig(savePath)
	#matplotlib.pyplot.show()			
			
def oldmain(): #lol!
	if argv[1].endswith(".txt"):
		bulkScaffoldPath = argv[1]
		scaffoldPaths = generateScaffoldFiles(argv[1])
	else:
		scaffoldPaths = []
		scaffoldPaths.append(os.getcwd() + "/" + argv[1])
	buildingBlocksPath = argv[2]
	linkersPath = argv[3]
	vdwRadPath = argv[4]
	smiLibPath = argv[5]

	#Creates each SMILE Library
	
	libraries = []
	for i in range(len(scaffoldPaths)):
		os.mkdir(scaffoldPaths[i][:(scaffoldPaths[i].rfind("/"))] + "/lib" + str(i))
		#print scaffoldPaths[i]
		SMILES.generateChemicals(scaffoldPaths[i], linkersPath, buildingBlocksPath, scaffoldPaths[i][:(scaffoldPaths[i].rfind("/"))] + "/lib" +str(i) +"/list" + str(i) + ".txt" , smiLibPath)
		#print "BUTTS"
		SMILES.isolateChemicals(scaffoldPaths[i][:(scaffoldPaths[i].rfind("/"))] + "/lib" +str(i) +"/list" + str(i) + ".txt", True, True, scaffoldPaths[i][:(scaffoldPaths[i].rfind("/"))] + "/lib" + str(i))
		libraries.append(scaffoldPaths[i][:(scaffoldPaths[i].rfind("/"))] + "/lib" + str(i))

	blastCap = "B"
	xOffset = 0
	yOffset = 0
	numSheests = 2
	graphs = []
	imageFiles = []
	textFiles = []
	#Generates all graphs
	for i in range(len(libraries)):
		##print libraries[i]
		files = sorted(os.listdir(libraries[i]))
		##print "HHHHHHHHHHHHHHHHHHHHHEEEEEEEY" + scaffoldPaths[i]
		os.mkdir(libraries[i] + "/baseline")
		baseline = generateBaseline(scaffoldPaths[i], vdwRadPath, blastCap, xOffset, yOffset, numSheests,libraries[i] + "/baseline/")
		##print files
		#print "Big Time go read"
		for j in range(len(files)):
			##print(files[j])
			if files[j].endswith(".xyz"):
				##print libraries[i] +"/" +files[j]
				
				currentGraph = subtractGraph(stacking.graphAnglebyOffset(libraries[i] + "/" +files[j], vdwRadPath, blastCap, xOffset, yOffset, numSheests), baseline)
				graphs.append(currentGraph)

				arrayGraphToHist(currentGraph, libraries[i] +"/" +  files[j][:-3] + ".png" )
				imageFiles.append(libraries[i] + "/" +files[j][:-3]  + ".svg ")
			
				textFiles.append(graphToText(currentGraph,libraries[i] + "/" +files[j][:-3]  + ".txt "))
				

def analyticalGraphFromFileToScatter(path, savePath, colorA):
	angle = []
	dist = []
	file = open(path, 'r')
	while True:
		try:
			readIn = file.next().split()
			#print readIn
		#	print readIn
			if readIn[2] > 0:
				angle.append(float(readIn[0]))
				dist.append(float(readIn[1]))
		except StopIteration:
			break
	file.close()
 	average = 0.0 
 	
 	for i in range(len(dist)):
 		average += dist[i]
 	average = average/len(dist)
 	SD = 0.0
 	for i in range(len(dist)):
 		SD += (dist[i] - average)*(dist[i] - average)
 	SD = SD/len(dist)
 	RSD = 100.0*(SD/average)
 	#print path
 	#print(len(angle))
 	#print(len(dist))
 	bminMaxList = getMinMaxList(angle,dist)
 	minMaxListA = bminMaxList[0]
 	minMaxListD = bminMaxList[1]
 	minMaxListI = bminMaxList[2]
 	baverage = 0.0 

 	
 	for i in range(len(minMaxListD)):
 		baverage += minMaxListD[i]
 	baverage = baverage/len(minMaxListD)
 	bSD = 0.0
 	for i in range(len(minMaxListD)):
 		bSD += (minMaxListD[i] - baverage)*(minMaxListD[i] - baverage)
 	bSD = bSD/len(minMaxListD)
 	bRSD = 100.0*(bSD/baverage)
 	#print  baverage, bSD, bRSD
 #	print baverage- bSD

 	minsA = []
 	minsD = []
 	#print minMaxListI
 	deepOneIndex = -1
 	deepTwoIndex = -1
 	toBeat = 200
 	for i in range(len(minMaxListD)):
 		if minMaxListD[i] < toBeat:
 			deepOneIndex = i
 			toBeat = minMaxListD[i]
 	toBeat = 200
 	for i in range(len(minMaxListD)):
 		if minMaxListD[i] < toBeat and i != deepOneIndex:
 			deepTwoIndex = i
 			toBeat = minMaxListD[i]

 	print minMaxListA[deepOneIndex],minMaxListA[deepTwoIndex]
 	if  0 - minMaxListA[deepOneIndex] >  0-minMaxListA[deepTwoIndex]:
 		furthestMin = minMaxListA[deepOneIndex]
 		closestMin = minMaxListA[deepTwoIndex]
 	else:
 		furthestMin = minMaxListA[deepTwoIndex]
 		closestMin = minMaxListA[deepOneIndex]
 	

 	stringA = 'Far Min: \n' + str(furthestMin) +" , " + str(dist[angle.index(furthestMin)]) + "\n " + str((dist[angle.index(furthestMin)]+1 +  dist[angle.index(furthestMin)]-1)/2)
 	stringB = 'Far Min: \n' + str(furthestMin) +" , " + str(dist[angle.index(furthestMin)]) + "\n " + str((dist[angle.index(furthestMin)]+1 +  dist[angle.index(furthestMin)]-1)/2)
 	stringC = 'Close Min: \n' + str(closestMin) +" , " + str(dist[angle.index(closestMin)]) + "\n " + str((dist[angle.index(closestMin)]+1 +  dist[angle.index(closestMin)]-1)/2)
	 
 	angle = []
	dist = []
	file = open(path, 'r')
	while True:
		try:
			readIn = file.next().split()
			#print readIn
		#	print readIn
			if readIn[2] > 0:
				angle.append(float(readIn[0]))
				dist.append(float(readIn[1]))
		except StopIteration:
			break
	file.close()


 	matplotlib.pyplot.scatter(angle,dist, color = colorA)
	matplotlib.pyplot.grid(b =True, which = 'major', linestyle = '-')
	matplotlib.pyplot.yticks(numpy.arange(3, 7+1, .5))
	matplotlib.pyplot.xticks(numpy.arange(-180, 180+1, 30))
	
	if furthestMin < 0:
		matplotlib.pyplot.annotate(stringA, xy=(furthestMin+30, 3.2), xycoords='data',bbox = dict(boxstyle="square,pad=0.3", fc ="white" ))
	
	else:
		matplotlib.pyplot.annotate(stringB, xy=(furthestMin-30, 3.2), xycoords='data',bbox = dict(boxstyle="square,pad=0.3", fc ="white" ))
	
	matplotlib.pyplot.annotate(stringC, xy=(closestMin, 3.2), xycoords='data',bbox = dict(boxstyle="square,pad=0.3", fc ="white" ))
	#matplotlib.pyplot.annotate('Average: ' + str(average) +" A\n" + "SD: " + str(SD) + " A \n" + "RSD: " + str(RSD) +" A", xy=(-180, 7), xycoords='data',bbox = dict(boxstyle="square,pad=0.3", fc ="white" ))
	#ax= matplotlib.pyplot.figure().add_subplot(111)
	#dist[angle.index(closestMin)]+1 +  dist[angle.index(closestMin)]-1)/2
	#ax.annotate('local max',  xytext=(2.5, -100), dict())


	#matplotlib.pyplot.show()
	matplotlib.pyplot.savefig(savePath)
	#matplotlib.pyplot.show()
	matplotlib.pyplot.clf()
	matplotlib.pyplot.cla()
	matplotlib.pyplot.close()

def digitalize(graph):
	for i in range(len(graph)):
		for j in range(len(graph[i])):
			if graph[i][j] > 0:
				graph[i][j] = 1
			else:
				graph[i][j] = 0

	return graph


def main():
	monomerPath = argv[1]
	linkersPath = argv[2]
	buildingBlocksPath = argv[3]
	smiLibPath = argv[4]
	vdwRadPath = argv[5]
	SMILES.grabAndSpin(argv[1])
	neomonomer = open(argv[1], 'r')

	#counter = 1
	#i = 2
	#currentGraph = stacking.graphAnglebyOffset(argv[1], vdwRadPath, 'B', 0, 0, 2)

	#arrayGraphToHist(currentGraph, "OUT.png")


















	#mainCentroid = stacking.centroid(neomonomer)
	ogLines = []
	woo = []
	waste = neomonomer.next()
	waste = neomonomer.next()
	os.mkdir("lib/")
	while True:
		try:
			#print "loop"
			readIn = neomonomer.next().split()
			ogLines.append(readIn)
		except StopIteration:
			#print "break"

			break
	counter = 0
	combolists = []
	for i in range(len(ogLines)):
		##print lines[i][0]
		lines = list(ogLines)
		if lines[i][0] == "H":
			os.mkdir("lib/site" + str(counter))
			#print "In der"
			lines[i][0] = "F"
			deathY = lines[i][1]
			hand = ""
			if lines[i][1] < 0:
				hand = "left"
			else:
				hand = "right"
			

			#print deathY
			for cc in range(len(lines)):
				
				#print float(lines[cc][2]) - float(deathY)
				
				if lines[cc][0] == "H" and float(lines[cc][2]) - float(deathY) < .001:
					
					lines = lines[:cc] + lines[cc+1:]
					break

			tempXYZ = open("lib/site" + str(counter)+ '/node' + str(counter) + "xyz.xyz", 'w')
			tempXYZ.write(str(len(lines))+ "\n")
			tempXYZ.write( "\n")

			for j in range(len(lines)):
				tempXYZ.write(str(lines[j][0]) + " " + str(lines[j][1]) + " "" " + str(lines[j][2]) + " "" " + str(lines[j][3]) + " \n") 
			tempXYZ.close()
			smile = SMILES.generateSMI("lib/site" + str(counter)+'/node' + str(counter) + "xyz.xyz")
			file = open(smile, 'r+')
			stringToEdit = ''
			while True:
				try:
					stringToEdit = str(file.next())
					##print stringToEdit
					stringToEdit = stringToEdit.replace("F", "([R1])")
					##print stringToEdit
				except StopIteration:
					break
			file.close()
			file = open(smile, 'w')
			file.write(stringToEdit.split()[0])
			file.close()
			outPath = "lib/site" + str(counter) + "/structureList" + str(counter) + ".smi"
			combolists.append(outPath)
			print "SMILING at Chemicals..."
			b=SMILES.generateChemicals(smile, linkersPath, buildingBlocksPath, outPath, smiLibPath)
			lines[i][0] = "H"
			print "Chasing those chemicals apart..."
			chems = SMILES.isolateChemicals(b,True,True,"lib/site" + str(counter),argv[1] ) 
			print "Here I am, rock you like a hurricane"
			print chems
			print hand
			for bb in range(len(chems)):
				print chems[bb][:-4] + ".xyz"
				SMILES.normalize(chems[bb][:-4] + ".xyz", stacking.centroid(argv[1]))
				bbcentroid = stacking.centroid(chems[bb][:-4] + ".xyz")
				print bbcentroid
				if hand == 'left':
					if bbcentroid[0] < .001:
						print "We good"
						waste = 0
					else:
						print "Spin Factor 1"
						SMILES.spinXY(chems[bb][:-4] + ".xyz", 180)
						#SMILES.generateSVGFromXYZ(chems[bb][:-4] + ".xyz")
				elif hand == 'right':
					if bbcentroid[0] > -.001:
						print "We good"
						waste = 0
					else:
						#SMILES.generateSVGFromXYZ(chems[bb][:-4] + ".xyz")
						print "Spin Factor 2"
						SMILES.spinXY(chems[bb][:-4] + ".xyz", 180)
				
				bbcentroid = stacking.centroid(chems[bb][:-4] + ".xyz")
				print bbcentroid
				SMILES.normalize(chems[bb][:-4] + ".xyz")
			blastCap = "B"
			xOffset = 0
			yOffset = 0
			numSheests = 2
			graphs= []
			cleanXYZTemp = open("lib/site" + str(counter)+'/node' + str(counter) + "xyz.xyz", "w")
			cleanXYZTemp.write(str(len(lines))+ " \n")
			cleanXYZTemp.write("$comment"+ " \n")
			for i in range(len(lines)):
				boo = lines[i]
				cleanXYZTemp.write(str(boo[0])+ " " + str(boo[1])+ " " + str(boo[2])+ " " + str(boo[3])+ " " + " \n")
			cleanXYZTemp.close()

			SMILES.grabAndSpin("lib/site" + str(counter)+'/node' + str(counter) + "xyz.xyz")
			smileYa = SMILES.generateSMI("lib/site" + str(counter)+'/node' + str(counter) + "xyz.xyz")
			print " Generating Baseline..."
			generateBaseline(monomerPath ,smileYa,vdwRadPath, blastCap, xOffset, yOffset, numSheests,"lib/site" + str(counter)+"/", argv[1])

			#baseline = digitalize( generateBaseline(monomerPath ,smileYa,vdwRadPath, blastCap, xOffset, yOffset, numSheests,"lib/site" + str(counter)+"/", argv[1]))
			#os.remove('UBERTEMP.xyz')
			
			for i in range(len(chems)):
				f = ""
				g = ""
				print "Graphing..."
				#currentGraph = subtractGraph(baseline, digitalize( stacking.graphAnglebyOffset(chems[i][:-3]+ "xyz", vdwRadPath, blastCap, xOffset, yOffset, numSheests)))
				f =stacking.precisionGraph(chems[i][:-4]+ ".xyz", vdwRadPath, .1, 0,0,2)
				g =coordsToText(f, "lib/site" + str(counter)+"/P.textGraph" + str(i) + ".txt")
				analyticalGraphFromFileToScatter("lib/site" + str(counter)+"/P.textGraph" + str(i) + ".txt", "lib/site" + str(counter)+"/P.graph" +str(i)+".png", "blue")

				
				#arrayGraphToHist(currentGraph, "lib/site" + str(counter)+"/graph" + str(i))
				"Histographing... :D"
				#graphToText(currentGraph, "lib/site" + str(counter)+"/textGraph" + str(i) + ".txt")		
				#c = graphFromFile("lib/site" + str(counter)+"/textGraph" + str(i) + ".txt")
				#print "lib/site" + str(counter)+"/textGraph" + str(i) + ".txt" + "SSSSSSSSSSSSSSSSss"
				#graphScatter(c, "lib/site" + str(counter)+"/graph" + str(i)+".png")
				
				#graphFromFileToScatter("lib/site" + str(counter)+"/textGraph" + str(i) + ".txt","lib/site" + str(counter)+"/graph" + str(i)+".png")
				#woo.append("lib/site" + str(counter)+"/textGraph" + str(i) + ".txt")

		
		counter += 1

'''	
base = digitalize(graphFromFile(argv[1]))
ss = digitalize(graphFromFile(argv[2]))

cc = subtractGraph(base, ss)


graphToText(cc, "SUB.txt")
graphScatter(cc, "SUB.png")
'''
main()











