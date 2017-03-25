#!/user/python/bin




from sys import argv
import os
import math
from scipy.spatial import distance
import sys
import numpy
import subprocess
import SMILES
import Centroid
class geomol:
	def __init__ (self, name, x,y,z):
		self.atmName = name
		self.atmX = x
		self.atmY = y
		self.atmZ = z
		###print(atmName)


#calculates the long axis on through the x-y plane, and recalculates the xy coordinates so that the y axis is the long axis

	
def centroid(filePath):
	'''
	toRead = open(filePath,'r',)
	totalAtoms = int(toRead.next())

	print toRead.next()
	xBar = 0.0;
	yBar = 0.0;
	zBar = 0.0;
	print range(totalAtoms)
	for line in toRead:
		lineList = line.split()
		xBar += float(lineList[1])
    	yBar += float(lineList[2])
    	zBar += float(lineList[3])
	
	for index in range (totalAtoms):
		lineList = toRead.next().split()
		
		print str(zBar) + "ZBAR1"
		xBar += float(lineList[1])
    	yBar += float(lineList[2])
    	zBar += float(lineList[3])
    	print str(zBar) + "ZBAR2"
    	##print lineList
    '''
	'''
	for index in range(0,2):
		###print index
		print index
		lineList = toRead.next().split()

    	xBar += float(lineList[1])
    	yBar += float(lineList[2])
    	zBar += float(lineList[3])
    	print "Line list"
    	print lineList
    	#print yBar
    	#print zBar
	
	xBar /= float(totalAtoms)
	yBar /= float(totalAtoms)
	zBar /= float(totalAtoms)

	array = []
	array.append(xBar)
	array.append(yBar)
	array.append(zBar)
	'''
	return Centroid.centroid(filePath)
	
def centroidFromArray(xArray, yArray, zArray):
	xBar = 0.0;
	yBar = 0.0;
	zBar = 0.0;
	totalAtoms = len(xArray)
	###print totalAtoms
	#print xArray
	for i in range(totalAtoms):
		#print xArray
		xBar += float(xArray[i])
    	yBar += float(yArray[i])
    	zBar += float(zArray[i])

	xBar /= float(totalAtoms)
	yBar /= float(totalAtoms)
	zBar /= float(totalAtoms)
	array = []
	array.append(xBar)
	array.append(yBar)
	array.append(zBar)
	return array
#Returns list of indeces of hydrogen atoms attached to acene carbons
def blast(blastCap, blastR, geomol):
	node = []
	
	aa = len(geomol.atmName)
	for i in range(aa):
		if(geomol.atmName[i] == blastCap):
			for j in range(len(geomol.atmName)):
				a = (geomol.atmX[i],geomol.atmY[i],geomol.atmZ[i])
				b = (geomol.atmX[j],geomol.atmY[j],geomol.atmZ[j])
				###print (a,b)
				###print(math.sqrt((int(geomol.atmX[j])-geomol.atmX[i])^2+(geomol.alignY[j]-geomol.alignY[i])^2+(geomol.alignZ[j]-geomol.alignZ[i])^2))
				#p#rint(isinstance(distance,int))
				###print(distance.euclidean(a,b))
				if distance.euclidean(a,b) > blastR and geomol.atmName[j] == 'H':
					###print("got")
					node.append(j)

	###print node
	return list(set(node))


def loadMoleculeTogeomol(monomerPath):
	newMol = reorienetMonomer(monomerPath)
	geo = geomol(newMol[0],newMol[1],newMol[2],newMol[3])


	return geo
	


#configFile = argv[1]#argv[1]
def buildFramework(monomerGeomol, blastCap, idealAngle, lateralSpace, xOffset, yOffset, numSheests):

	#loads calculations parameters from config file
	#	Format
	# 		Path to monomer to be stacked
	#       Ideal angle between monomer
	#		Blast cap elements
	#		(Optional) Pattern: "zigzag", "insquare", "outsquare"

	
	

	framework = []
	newatm = monomerGeomol
	#print newatm[0],newatm[1],newatm[2],newatm[3]
	framework.append(newatm)

	for i in range(1,numSheests):
		###print(i)
		newName = []
		newXa = []
		newYa = []
		newZa = []
		#print framework[i-1].atmX
		centroidData =  centroidFromArray(framework[i-1].atmX, framework[i-1].atmY, framework[i-1].atmZ)
		centroidX = centroidData[0]
		centroidY = centroidData[1]
		centroidZ = centroidData[2]
		xCorrect = framework[0].atmX[0]
		zCorrect = framework[0].atmZ[0]
		#print xCorrect
		#print zCorrect
		'''
		for i1 in range(len(framework)):
			for j1 in range(len(framework[i-1].atmName)):
				framework[i1].atmX[j1]  -= centroidX + xCorrect
				framework[i1].atmY[j1]  -= centroidY
				framework[i1].atmZ[j1] -= centroidZ + zCorrect
				
		
		'''


		for j in range(len(framework[i-1].atmName)):

			
			newX = (framework[i-1].atmX[j]*numpy.cos(numpy.deg2rad(idealAngle)) + framework[i-1].atmZ[j]*numpy.sin(numpy.deg2rad(idealAngle))) + xOffset
  #(framework[i-1].atmX[j]*numpy.cos(numpy.deg2rad(idealAngle)) - framework[i-1].atmZ[j]*numpy.sin(numpy.deg2rad(idealAngle))) + xOffset  #(framework[i-1].atmX[j]*numpy.cos(numpy.deg2rad(idealAngle)) - framework[i-1].atmZ[j]*numpy.sin(numpy.deg2rad(idealAngle))) + xOffset 
			#framework[i-1].atmX[j]+ xOffset

			#(framework[i-1].atmX[j]*numpy.cos(numpy.deg2rad(idealAngle)) - framework[i-1].atmZ[j]*numpy.sin(numpy.deg2rad(idealAngle))) + xOffset
	
			newY = framework[i-1].atmY[j]+ yOffset
			 #framework[i-1].atmY[j] + yOffset

			newZ = (framework[i-1].atmZ[j]*numpy.cos(numpy.deg2rad(idealAngle)) - framework[i-1].atmX[j]*numpy.sin(numpy.deg2rad(idealAngle))) + lateralSpace # framework[i-1].atmZ[j]+ lateralSpace# (framework[i-1].atmZ[j]*numpy.cos(numpy.deg2rad(idealAngle)) + framework[i-1].atmX[j]*numpy.sin(numpy.deg2rad(idealAngle))) + lateralSpace
			
			
			newName.append(framework[i-1].atmName[j])
			newXa.append(float("{:10.4f}".format(newX)))
			newYa.append(float("{:10.4f}".format(newY)))
			newZa.append(float("{:10.4f}".format(newZ)))
		newatm = geomol(newName,newXa, newYa,newZa)
		framework.append(newatm)

	return framework



class vdwMol:
	def __init__ (self, name, x,y,z,rad):
		self.atmName = name
		self.atmX = x
		self.atmY = y
		self.atmZ = z
		self.atmRad = rad
	

def loadvdwData(vdwDataPath):
	atmName = []
	atmRadii = []
	readIn = open(vdwDataPath, 'r')

	while True:

		try:
			line = readIn.next().split()
			atmName.append(line[0])
			atmRadii.append(line[1])



		except StopIteration:
			break



	returnData = []
	returnData.append(atmName)
	returnData.append(atmRadii)
	return returnData
def frameworkTovdw(framework):
	vdwFramework = []
	###print len(framework)
	name = ['C', 'H', "O", 'B', 'N', 'S']
	data = [1.7,1.2,1.5,1.9,1.6,1.8]
	for i in range (len(framework)):
		
		sheetvdwRadii = []
		newMol = ""
		for j in range(len(framework[i].atmName)):
			
			for k in range(len(name)):
				if framework[i].atmName[j] == name[k]:
					sheetvdwRadii.append(data[k])
					#print "Serious print"
					#print framework[i].atmName[j]
					#print data[k]
					break
		newMol = vdwMol(
			framework[i].atmName,
			framework[i].atmX,
			framework[i].atmY,
			framework[i].atmZ,
			sheetvdwRadii)
		
		vdwFramework.append(newMol)
	###print len(vdwFramework)
	return vdwFramework


def distanceMeasureSingleValue(framework):
	###print len(framework)
	minDist = "init"
	minIndexPair = []
	for i in range(len(framework[0].atmName)):
		for j in range(len(framework[1].atmName)):
				a = (framework[0].atmX[i],framework[0].atmY[i],framework[0].atmZ[i])
				b = (framework[1].atmX[j],framework[1].atmY[j],framework[1].atmZ[j])
				dist =  distance.euclidean(a,b) - framework[0].atmRad[j] - framework[1].atmRad[j]
				###print dist
				if minDist == "init":
					###print dist
					minDist = dist
					minIndexPair = [i,j]
				elif  dist < minDist:
					###print "|   " + str(dist)
					minDist = dist
					minIndexPair = [i,j]


	
		##print "Closest Atoms:"
	##print minDist, minIndexPair

	returnData = []
	returnData.append(minDist)
	returnData.append(minIndexPair)
	return returnData

def distanceMeasureAllValue(framework):
	###print len(framework)
	collisionDist = []
	collisionIndexPair = []
	for i in range(len(framework[0].atmName)):
		for j in range(len(framework[1].atmName)):
				a = (framework[0].atmX[i],framework[0].atmY[i],framework[0].atmZ[i])
				b = (framework[1].atmX[j],framework[1].atmY[j],framework[1].atmZ[j])
				dist =  distance.euclidean(a,b) - framework[0].atmRad[i] -  framework[0].atmRad[j]
				#print dist
				if dist < 0:
					collisionDist.append(dist)
					collisionIndexPair.append([i,j])
					#print "Collision#: " + str(len(collisionIndexPair)) +" Line: " + str(i+1),  framework[0].atmName[i], framework[0].atmRad[i], "Line: " + str(len(framework[0].atmName) + j+1), framework[1].atmName[j],  framework[0].atmRad[j]
					#print framework[0].atmX[i],framework[0].atmY[i],framework[0].atmZ[i]
					#print framework[1].atmX[j],framework[1].atmY[j],framework[1].atmZ[j]
					#print distance.euclidean(a,b)
					#print dist
					#print  distance.euclidean(a,b) - framework[0].atmRad[i] -  framework[0].atmRad[j]
					
	returnData = []
	returnData.append(collisionDist)
	returnData.append(collisionIndexPair)
	
	#if len(collisionDist) < 1:
	#	print "No collision"
	#	print collisionDist
	#else:
	#	print "		" + str(len(collisionDist)) + " Collisions"
	#	print collisionDist
		
	return returnData



def writeXYZ(framework, savePath):
	testxyz = open( savePath, 'w')
	testxyz.write(str(len(framework)*(len(framework[0].atmName))) + '\n' + '\n')
	for i in range (len(framework)):
		for j in range (len(framework[i].atmName)):
			#if j == 0 :
				###print framework[i].atmName[j] + " " + str(framework[i].atmX[j]) + " " + str(framework[i].atmY[j]) + " " + str(framework[i].atmZ[j])
			testxyz.write(str(framework[i].atmName[j]) + " " + str(framework[i].atmX[j]) + " " + str(framework[i].atmY[j]) + " " + str(framework[i].atmZ[j]) + "\n")
	testxyz.close()



def precisionGraph(monomerPath, vdwDataPath, errorMargin, xOffset,yOffset,numSheests):
	coords = []
	print "		Precision Graph enter: "
	for iAngle in range(0,360):
		coords.append([iAngle,searchContact(monomerPath, vdwDataPath, iAngle, errorMargin, xOffset,yOffset,numSheests)])
	#print coords
	print "		Precision graph exit: "
	return coords
def searchContact(monomerPath, vdwDataPath, angle, errorMargin, xOffset,yOffset,numSheests):
	
	ss =  open(monomerPath, "r")
	ssName = []
	ssX = []
	ssY = []
	ssZ = []
	waste = ss.next()
	waste = ss.next()
	while True:
		try:
			readIn = ss.next().split()
			ssName.append(readIn[0])
			ssX.append(float(readIn[1]))
			ssY.append(float(readIn[2]))
			ssZ.append(float(readIn[3]))
		except StopIteration:
			break
	z = geomol(ssName, ssX, ssY, ssZ)
	dist = 0
	farthestContact = 0
	while True:
		a = buildFramework(z, "", angle, dist, xOffset, yOffset, numSheests)
		c = frameworkTovdw(a)

		d = distanceMeasureAllValue(c)
		if len(d[0]) > 0:
			farthestContact = dist
			dist += .01
		
		else:
			break
		#print farthestContact, dist
	if farthestContact > 0:
		#print "RUN"
		return recursiveNarrow(monomerPath, vdwDataPath, angle, errorMargin, xOffset,yOffset,numSheests, float(farthestContact), float(dist))
	else:
		#print "NOOOO"
		return 0
def recursiveNarrow(monomerPath, vdwDataPath, angle, errorMargin, xOffset,yOffset,numSheests, contact, hole):
	#print hole, contact
	#print contact + ((hole - contact)/2)
	#print "----------------"
	ss =  open(monomerPath, "r")
	ssName = []
	ssX = []
	ssY = []
	ssZ = []
	waste = ss.next()
	waste = ss.next()
	while True:
		try:
			readIn = ss.next().split()
			ssName.append(readIn[0])
			ssX.append(float(readIn[1]))
			ssY.append(float(readIn[2]))
			ssZ.append(float(readIn[3]))
		except StopIteration:
			break
	z = geomol(ssName, ssX, ssY, ssZ)
	a = buildFramework(z, "", angle, hole - (hole - contact)/2, xOffset, yOffset, numSheests)
	c = frameworkTovdw(a)
	d = distanceMeasureAllValue(c)
	if len(d[0]) > 0:
		#print hole - (hole - contact)/2
		if hole - contact + (hole - contact)/2.0 <= errorMargin:
			#print "HERE" #+ str((hole - contact)/2.0)
			return contact +  (hole - contact)/2.0
		else:
			#print "THERE" #+ str((hole - contact)/2.0)
			return recursiveNarrow(monomerPath, vdwDataPath, angle, errorMargin, xOffset,yOffset,numSheests,contact + (hole - contact)/2.0,hole)
	else:
		#print "WHERE" #+ str((hole - contact)/2.0)
		return recursiveNarrow(monomerPath, vdwDataPath, angle, errorMargin, xOffset,yOffset,numSheests,contact, contact + (hole - contact)/2.0)
def graphAnglebyOffset(monomerPath, vdwDataPath,  blastCap, xOffset, yOffset, numSheests, deposit = False, savePath = "NA"):
	##print "Entering graphAnglebyOffset"
	graph = []
	print "		You spin me round!"

	ss =  open(monomerPath, "r")
	ssName = []
	ssX = []
	ssY = []
	ssZ = []
	waste = ss.next()
	waste = ss.next()
	while True:
		try:
			readIn = ss.next().split()
			ssName.append(readIn[0])
			ssX.append(float(readIn[1]))
			ssY.append(float(readIn[2]))
			ssZ.append(float(readIn[3]))
		except StopIteration:
			break
	z = geomol(ssName, ssX, ssY, ssZ)
	print "		Making a bunch of frameworks..."
	for iAngle in range(0,360):
	#	iAngle = 258
		row = []
		''
		for jDistance in numpy.arange(4,10,0.25):
			
			#jDistance = 6.0
			#print "					DISTANCE    :    "  + str(jDistance)
			a = buildFramework(z, blastCap, iAngle, jDistance, xOffset, yOffset, numSheests)
			#if deposit and iAngle > 240 and iAngle < 325 and jDistance > 5.75 and jDistance < 6.5:
			#	writeXYZ(a,savePath + "/" + monomerPath[:-3] + str(iAngle) +"-" + str(jDistance)+".xyz")

			#print "HEEE"
			#writeXYZ(a,"AAA" + str(iAngle) + " " + str(jDistance) + ".xyz")
			#b = loadvdwData(vdwDataPath)
			#print "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE     : " 
			#print b
			c = frameworkTovdw(a)#3, b[0], b[1])

			d = distanceMeasureAllValue(c)
			#if iAngle == 225 and jDistance == 7.75:
				#print len(d[0])
			#point = 0
			#print len(d[0])
			#if len(d[0]) > 0:
				#point = 1
			point =  len(d[0]) 

			row.append(point)
			#print row
		#print row
		graph.append(row)
		###print str(iAngle) + " " + str(column)
	##print "bye"
	return graph



def main():
#monomerPath, blastCap, idealAngle, lateralSpace, xOffset, yOffset, numSheests
	writeXYZ(buildFramework(argv[1], 'B', 0, 3, 0,0,1))

#main()


