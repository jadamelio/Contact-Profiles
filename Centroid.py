#!/usr/bin/python

from sys import argv


def centroid(Ar):
	toRead = open(Ar,'r')

	totalAtoms = int(toRead.next())

	toRead.next()
	xBar = 0.0;
	yBar = 0.0;
	zBar = 0.0;
	for index in range(totalAtoms):
	    lineList = toRead.next().split()
	    xBar += float(lineList[1]);
	    yBar += float(lineList[2]);
	    zBar += float(lineList[3]);

	xBar /= totalAtoms;
	yBar /= totalAtoms;
	zBar /= totalAtoms;

	#print xBar, yBar, zBar

	array = []
	array.append(xBar)
	array.append(yBar)
	array.append(zBar)
	return array
	toRead.close()
	#return 0


#centroid(argv[1])





