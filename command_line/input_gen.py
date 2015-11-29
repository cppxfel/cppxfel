distance = 0
wavelength = 0
centre = (0, 0)
pixelSize = (0, 0)
spacegroup = 0
unit_cell_dimensions = (0, 0, 0, 0, 0, 0)
panels = []

"""
Function finds the filename stem from an experiment JSON file
excluding the extension.

Args: filename: string filename of JSON file.

Returns: filename stem.

"""
def rootnameStem(filename):
	return filename.replace("_experiments.json", "").replace("_strong.list", "")

"""
Function prints matrices from DIALS experiment objects separating
into unitcell and rotation components into a string stream for a
single experiments JSON file.

Args:	experiments: DIALS experiment object list
			filename: Name of experiments JSON file (also source of
			name of image)
			output: string IO object

Returns: None

"""
def printExperiments(experiments, filename, output):	
	rootname_stem = rootnameStem(filename)
	print >>output, "image " + rootname_stem[1:]
	
	global distance, centre, wavelength, pixelSize, spacegroup, unit_cell_dimensions
	
	if len(experiments):
		beam = experiments[0].beam
		detector = experiments[0].detector[0]
		distance = detector.get_distance()
		centre = detector.get_ray_intersection_px(beam.get_s0())
		wavelength = beam.get_wavelength()

		pixelSize = detector.get_pixel_size()

		print >> output, "wavelength", wavelength
		print >> output, "distance", distance
		print >> output, "centre", centre[0], centre[1]

	spotListName = rootname_stem + "_strong.list"
	
	if (os.path.isfile(spotListName)):
		print >> output, "spots", spotListName

	for experiment in experiments:
		spacegroup = experiment.crystal.get_space_group().type().number()
		unit_cell_dimensions = experiment.crystal.get_unit_cell().parameters()
		unit_cell = experiment.crystal.get_B()	
		rotation = experiment.crystal.get_U()	

		print >>output, "unitcell",

		for component in unit_cell:
			print >>output, component,

		print >> output, ""

		print >>output, "rotation",

		for component in rotation:
			print >>output, component,

		print >> output, ""


"""
Takes an individual experiments JSON file and runs DIALS
to retrieve the matrix objects, and then sends the matrices off
to be printed to the string IO.

args:	filename: name of experiments JSON file
			output: string IO object

returns: None

"""
def matrixForFilename(filename, output):
	arguments = []
	arguments.append(filename)

	from dials.util.options import OptionParser
	from dials.util.options import flatten_experiments
	from libtbx.utils import Abort

	parser = OptionParser(
	read_experiments=True,
	check_format=False)

	params, options = parser.parse_args(args=arguments, show_diff_phil=True)
	experiments = flatten_experiments(params.input.experiments)

	path = experiments[0].imageset.get_path(0)
	print "Image", path, "has", len(experiments), "experiments."
	paths.append(path)
	printExperiments(experiments, filename, output)

def printSpots():
	import glob
	global output
	globString = "*_strong.list"
	files = glob.glob(globString)
	
	for file in files:
		rootname = rootnameStem(file)
		
		if os.path.isfile(rootname + "_experiments.json"):
			continue
		
		imagePath = rootname[1:] + ".pickle"
		
		paths.append(imagePath)
		printExperiments([], file, output)

def printMatrices():
	import glob
	global output
	globString = "*_experiments.json"
	files = glob.glob(globString)
	
	for file in files:
		matrixForFilename(file, output)	

def dumpPanels(image):
	global panels
	print "Dumping panels from first image"
	f = open(image, 'rb')
	x = pickle.load(f)
	panels = x['ACTIVE_AREAS']
	panelTxt = StringIO.StringIO()
	
	for i in range(0, len(panels), 4):
		print >> panelTxt, "PANEL",
		print >>panelTxt, panels[i + 1], panels[i + 0], panels[i + 3], panels[i + 2]
	
	outputFilename = "panels.txt"
	file = open(outputFilename, 'w')
	print >>file, panelTxt.getvalue()
	file.close()

def dumpImages(imagePaths):
	for path in imagePaths:
		print "Dumping image:", path
		rootname = splitext(basename(path))[0]
		db = dxtbx.load(path).get_detectorbase()
		
		data = db.get_raw_data()

		data_array = array.array('i')

		for i in range(0, len(data)):
			data_array.append(data[i])

		string = data_array.tostring()
		newFile = open(rootname + '.img', 'wb')
		newFile.write(string)
		newFile.close()

from os.path import basename
from os.path import splitext
import StringIO
import dxtbx
import pickle
import array
import os
import sys
from multiprocessing import Process
from scitbx.array_family import flex
import scitbx.math

paths = []

output = StringIO.StringIO()

printMatrices()
printSpots()

print output.getvalue()

outputFilename = "matrices.dat"
outputFile = open(outputFilename, 'w')
print >>outputFile, output.getvalue()[:-1]
outputFile.close()

threads = []
thread_count = int(os.getenv('NSLOTS', 4))
image_num = len(paths)
images_per_thread = image_num / thread_count

if (len(paths) == 0):
	print "No images with matrices provided."
	exit()

print "Dumping panels from first image..."
dumpPanels(paths[0])

print "Total images ready for dumping: ", image_num

for i in range(0, thread_count):
	min = int(i * images_per_thread)
	max = int((i + 1) * images_per_thread)

	print "Dumping", min, "to", max, " images on this thread"
	images = paths[min:max]
	
	thread = Process(target=dumpImages, args=(images, ))
	threads.append(thread)
	thread.start()

print "Creating input files"

integrateTxt = StringIO.StringIO()

print >> integrateTxt, "ORIENTATION_MATRIX_LIST matrices.dat"
print >> integrateTxt, "NEW_MATRIX_LIST orientations.dat\n"

if spacegroup > 0:
	print >> integrateTxt, "SPACE_GROUP", spacegroup

if len(unit_cell_dimensions):
	print >> integrateTxt, "UNIT_CELL",
	for dimension in unit_cell_dimensions:
		print >> integrateTxt, dimension,
print >> integrateTxt 
print >> integrateTxt, "FIX_UNIT_CELL ON\n"
print >> integrateTxt, "MM_PER_PIXEL", pixelSize[0]
print >> integrateTxt, "BEAM_CENTRE", centre[0], centre[1]
print >> integrateTxt, "DETECTOR_DISTANCE", distance
print >> integrateTxt, "INTEGRATION_WAVELENGTH", wavelength, "\n"
print >> integrateTxt, "PANEL_LIST panels.txt"
print >> integrateTxt, "METROLOGY_SEARCH_SIZE 1\n"
print >> integrateTxt, "SHOEBOX_FOREGROUND_PADDING 1"
print >> integrateTxt, "SHOEBOX_NEITHER_PADDING 2"
print >> integrateTxt, "SHOEBOX_BACKGROUND_PADDING 3\n"
print >> integrateTxt, "INTENSITY_THRESHOLD 12"
print >> integrateTxt, "ABSOLUTE_INTENSITY OFF"
print >> integrateTxt, "REFINE_ORIENTATIONS ON\n"
print >> integrateTxt, "VERBOSITY_LEVEL 1\n"
print >> integrateTxt, "COMMANDS\n"
print >> integrateTxt, "INTEGRATE"

outputFilename = "integrate.txt"
outputFile = open(outputFilename, 'w')
print >>outputFile, integrateTxt.getvalue()
outputFile.close()

print "New template input file integrate.txt"
print "Please check your target unit cell and space group"

refineTxt = StringIO.StringIO()

print >> refineTxt, "ORIENTATION_MATRIX_LIST refine-orientations.dat"
print >> refineTxt, "MATRIX_LIST_VERSION 2.0"
print >> refineTxt, "NEW_MATRIX_LIST refined.dat\n"
print >> refineTxt, "\nPARTIALITY_CUTOFF 0.2\n"
print >> refineTxt, "COMMANDS\n"
print >> refineTxt, "REFINE_PARTIALITY"

outputFilename = "refine.txt"
outputFile = open(outputFilename, 'w')
print >>outputFile, refineTxt.getvalue()
outputFile.close()

print "New template input file refine.txt"

mergeTxt = StringIO.StringIO()

print >> mergeTxt, "ORIENTATION_MATRIX_LIST merge-orientations.dat"
print >> mergeTxt, "MATRIX_LIST_VERSION 2.0"
print >> mergeTxt, "INITIAL_MTZ allMerge5.mtz\n"
print >> mergeTxt, "OUTLIER_REJECTION_SIGMA 1.6"
print >> mergeTxt, "MERGE_ANOMALOUS OFF"
print >> mergeTxt, "RECALCULATE_SIGMA ON"
print >> mergeTxt, "CORRELATION_THRESHOLD 0.9"
print >> mergeTxt, "COMMANDS\n"
print >> mergeTxt, "MERGE"

outputFilename = "merge.txt"
outputFile = open(outputFilename, 'w')
print >>outputFile, mergeTxt.getvalue()
outputFile.close()

print "New template input file merge.txt"

if distance == 0:
	print "You MUST change the experimental parameters in integrate.txt"
