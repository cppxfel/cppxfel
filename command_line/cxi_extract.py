#!/usr/bin/env libtbx.python

import h5py
import pickle
import sys
import scitbx_array_family_flex_ext
import numpy
import os
from multiprocessing import Process
from scitbx.array_family.flex import grid


filename = ""
geometry = ""
entry = "entry_1"
beamX = 882.
beamY = 882.
width = 1765
height = 1765
skip = 0

print "CXI to pickle dumper"

for arg in sys.argv[1:]:
	key_value = arg.split("=")
	if (len(key_value) < 2):
		print "Did not understand", arg
		continue
	
	if (key_value[0] == "filename"):
		filename = key_value[1]
	elif (key_value[0] == "entry"):
		entry = key_value[1]
	elif (key_value[0] == "geometry"):
		geometry = key_value[1]
	elif (key_value[0] == "beamx"):
		beamX = float(key_value[1])
	elif (key_value[0] == "beamy"):
		beamY = float(key_value[1])
	elif (key_value[0] == "width"):
		width = int(key_value[1])
	elif (key_value[0] == "height"):
		height = int(key_value[1])
	elif (key_value[0] == "skip"):
		skip = int(key_value[1])
	else:
		print "Did not understand", arg

if len(filename) == 0:
	print "No filename selected, use filename=xxx.cxi"
	exit()
	
if len(geometry) == 0:
	print "No geometry file selected, use geometry=xxx.geom"
	exit()
	
print "Selected filename:", filename
print "Dumping entry:", entry
print "Applying geometry:", geometry

cxi = h5py.File(filename, 'r')

num_images = cxi[entry + "/data_1/data"].len()

print "Examining geometry file: finding rigid_group_d0."
panels = []
geomfile = open(geometry, 'r')
beginning = geomfile.tell()

for line in geomfile:
	if line[0] == ";":
		continue
	components = line.split("=")
	if (components < 2):
		continue
	if components[0].strip() == "rigid_group_d0":
		print "Found rigid_group_d0"
		panels = components[1].strip().split(",")
		print "This group has", len(panels), "panels."
		break

panelInfo = {}

for panel in panels:
	geomfile.seek(beginning)
	panelInfo[panel] = {}
	for line in geomfile:
		if (line[:len(panel)] == panel and line[len(panel):len(panel)+1] == "/"):
			key_value = line.split("=")
			key = key_value[0].strip()
			value = key_value[1].strip()
			
			id_property = key.split("/")
			id = id_property[0]
			property = id_property[1]
			
			if len(value.split(" ")) > 1:
				subProperties = value.split(" ")
				for i in range(len(subProperties)):
					subProperties[i] = subProperties[i].replace("x", "").replace("y", "")
					subProperties[i] = float(subProperties[i])
				value = subProperties
			else:
				value = float(value)
			
			print "Found value", value, "for panel", panel, ", property", property
			panelInfo[panel][property] = value

geomfile.close()

print "\nEntry has", num_images, "images.\n"
sample_file = open('sample.pickle', 'rb')
sample = pickle.load(sample_file)

for i in range(skip, num_images):
	image = cxi[entry + "/data_1/data"][i]
	identifier = cxi[entry + "/data_1/experiment_identifier"][i]
	distance = cxi[entry + "/data_1/distance"][i] * 1000
	pixelSize = cxi[entry + "/data_1/x_pixel_size"][i] * 1000
	energy = cxi[entry + "/instrument_1/source_1/energy"][i]
	keV = energy * 6.2416006565e+15
	wavelength = 12.4 / keV
	print "Dumping", identifier
	alldata = []
	count = 0
	rowSize = len(image[0])
	for row in image:
		rowdata = numpy.array(row, dtype=numpy.int32)
		alldata.extend(row)
		count += 1
	totalPicklePixels = width * height
	print "Generating canvas for pickle image with", totalPicklePixels, "pixels."

	picklePixels = scitbx_array_family_flex_ext.int()
	picklePixels.resize(totalPicklePixels)

	for panel in panelInfo:
		origTopLeftX = int(panelInfo[panel]['min_fs'])
		origTopLeftY = int(panelInfo[panel]['min_ss'])
		origBottomRightX = int(panelInfo[panel]['max_fs'])
		origBottomRightY = int(panelInfo[panel]['max_ss'])
		relativeNewX = panelInfo[panel]['corner_x']
		relativeNewY = panelInfo[panel]['corner_y']
		absoluteNewX = int(beamX + relativeNewX)
		absoluteNewY = int(beamY + relativeNewY)
		axisXX = int(round(panelInfo[panel]['fs'][0]))
		axisXY = int(round(panelInfo[panel]['fs'][1]))
		axisYX = int(round(panelInfo[panel]['ss'][0]))
		axisYY = int(round(panelInfo[panel]['ss'][1]))

		print panel + "\t" + str(origTopLeftX) + "\t" + str(origTopLeftY) + "\t" + str(origBottomRightX) + "\t" + str(origBottomRightY) + "\t" + str(absoluteNewX) + "\t"  + str(absoluteNewY) + "\t" + str(axisXX) + "\t" + str(axisXY) + "\t" + str(axisYX) + "\t" + str(axisYY)

		for k in range(origTopLeftX, origBottomRightX):
			xOffset = k - origTopLeftX
			for j in range(origTopLeftY, origBottomRightY):
				yOffset = j - origTopLeftY
				newX = absoluteNewY + axisXX * yOffset + axisXY * xOffset
				newY = absoluteNewX + axisYX * yOffset + axisYY * xOffset
			
				newValue = int(alldata[j * rowSize + k])
				picklePixels[newX * width + newY] = newValue
	
	string = picklePixels.as_numpy_array().tostring()
	newFile = open(identifier + '.img', 'wb')
	newFile.write(string)
	newFile.close()

	sample['DISTANCE'] = distance
	sample['SIZE1'] = width
	sample['SIZE2'] = height
	sample['ACTIVE_AREAS'] = scitbx_array_family_flex_ext.int([0, 0, width, height])
	sample['PIXEL_SIZE'] = pixelSize
	sample['BEAM_CENTER_X'] = 97.02
	sample['BEAM_CENTER_Y'] = 97.02
	sample['WAVELENGTH'] = wavelength

	picklePixels.reshape(grid(width, height))

	sample['DATA'] = picklePixels

	pickleName = identifier + '.pickle'
	new_pickle = open(pickleName, 'wb')
	pickle.dump(sample, new_pickle)
	print "Dumped pickle", pickleName, "and img", identifier + ".img"

	newFile.close()
	new_pickle.close()


print "Done."