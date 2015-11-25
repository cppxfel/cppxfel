from __future__ import division

import boost.python
from os import path
ext = boost.python.import_ext("cppxfel_ext")
from cppxfel_ext import *

class Parser:

	def __init__(self, commandFile):
		self.parser = ext.cppParser(commandFile)

# file provided under keyword ORIENTATION_MATRIX_LIST
	def loadImagesFromFile(self):
		

# 	def loadImages(self, experimentList):
# 		import dxtbx
# 		import array
# 	
# 		firstExp = experimentList[0]
# 		imageNameLong = str(firstExp.imageset)[2:-2]
# 		imageName = path.basename(imageNameLong)
# 		db = dxtbx.load(imageName).get_detectorbase()
# 		imageData = db.get_raw_data()
# 		data_array = array.array('i')
# 		
# 		for i in range(0, len(imageData)):
# 			data_array.append(imageData[i])
# 			
# 		length = len(imageData)
# 		distance = firstExp.detector[0].get_distance()
# 		wavelength = firstExp.beam.get_wavelength()
# 		self.parser.loadImage(imageName, data_array.tostring(), distance, wavelength)
# 			
# 		for experiment in experimentList:
# 			unit_cell = list(experiment.crystal.get_B().as_mat3())
# 			rotation = list(experiment.crystal.get_U().as_mat3())
# 			self.parser.addMatrixToLastImage(unit_cell, rotation)

	def integrate(self):
		self.parser.integrate()
		