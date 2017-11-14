#!/usr/bin/env cppxfel.python

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

def baseName(filename):
    return os.path.splitext(filename)[0]

"""
Function prints matrices from DIALS experiment objects separating
into unitcell and rotation components into a string stream for a
single experiments JSON file.

Args:   experiments: DIALS experiment object list
                        filename: Name of experiments JSON file (also source of
                        name of image)
                        output: string IO object

Returns: None

"""
def printExperiments(experiments, filename, output):
    rootname_stem = rootnameStem(filename)

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

args:   filename: name of experiments JSON file
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
    printExperiments(experiments, filename, output)


"""
Goes through each image path provided and searches for an associated
_*strong.list file and _*experiments.json file. Prints relevant data
if the associated file is found.
arg: path to an individual file.

"""
def printData(path):
    global output
    baseFile = baseName(path)
    spotPickleName = "_" + baseFile + "_strong.pickle"
    spotListName = spotPickleName.replace('pickle', 'list')
    experimentsJson = "_" + baseFile + "_experiments.json"
    from cppxfel import dials_print_spots

    print "Finding data for", path

    if (os.path.isfile(path)):
        print >> output, "image", baseFile

    if os.path.isfile(spotPickleName):
        print "Found spots."
        dials_print_spots.printSpots([spotPickleName])

    if (os.path.isfile(spotListName)):
        print >> output, "spots", spotListName
    else:
        print >> output, "spots find"

    if os.path.isfile(experimentsJson):
        print "Found DIALS indexing solutions."
        matrixForFilename(experimentsJson, output)

def dumpPanels(images):
    import pickle
    global panels
    global distance, wavelength, centre, pixelSize
    print "Attempting to dump panel file"
    # i is an integer specifying the image counter
    i = 0
    successful = False

    while not successful:
        try:
            f = open(images[i], 'rb')
            x = pickle.load(f)
        except Exception as e:
            print Exception
            i += 1

            if (i > len(images)):
                break

            continue

        # now we dump from this pickle file and stop the while loop

        successful = True

        try:
            panels = x['ACTIVE_AREAS']
            distance = x['DISTANCE']
            wavelength = x['WAVELENGTH']
            pixelSize = (x['PIXEL_SIZE'], 0)
            centre = (x['BEAM_CENTER_X'] / pixelSize[0], x['BEAM_CENTER_Y'] / pixelSize[0])
        except Exception as e:
            print Exception
            break

        # TODO: find distance, wavelength, beam centre from first file.

        panelTxt = StringIO.StringIO()

        for i in range(0, len(panels), 4):
            print >> panelTxt, "PANEL",
            print >>panelTxt, panels[i + 1], panels[i + 0], panels[i + 3], panels[i + 2]

        outputFilename = "panels.txt"
        file = open(outputFilename, 'w')
        print >>file, panelTxt.getvalue()
        file.close()
        return

    # If we got here, a pickle file was not found.

    print "No pickle files found. Please check your panels.txt file."
    outputFilename = "panels.txt"
    file = open(outputFilename, 'w')
    print >>file, "# Enable for basic CSPAD detector."
    print >>file, "#PANEL 0 0 1765 1765"
    print >>file, "# Enable for Rayonix detector."
    print >>file, "#PANEL 0 0 1920 1920"
    file.close()
    return

"""
        Dumps images specified in the imagePath list from
        pickle form into img form, but will skip those
        which have already been dumped.
        arg imagePaths: list of strings.
"""
def dumpImages(imagePaths):
    for path in imagePaths:
        print "Dumping image:", path
        rootname = splitext(basename(path))[0]
        imageName = rootname + '.img'

        if (os.path.exists(imageName)):
            print "Image", imageName, "already dumped - skipping."
            continue

        ext = splitext(basename(path))[1]

        if not ext == ".pickle":
            print "Warning: not a pickle file. Converting to pickle first."
            command = "cxi.image2pickle " + path
            os.system(command)

        db = dxtbx.load(rootname + ".pickle").get_detectorbase()
        data = db.get_raw_data()

        data_array = array.array('i')

        for i in range(0, len(data)):
            data_array.append(data[i])

        string = data_array.tostring()
        newFile = open(imageName, 'wb')
        newFile.write(string)
        newFile.close()

# **********************************
# ***** START OF EXECUTED CODE *****
# **********************************

# Take all the paths provided after cppxfel.input_gen.
paths = sys.argv[1:]

if (len(paths) == 0):
    print "Please specify the images you want to dump, e.g.:"
    print "\tcppxfel.input_gen shot*.pickle"
#       print "\tcppxfel.input_gen shot*.cbf"
    exit()

# Start the output stringIO which will be written to matrices.dat
output = StringIO.StringIO()

# For each image, print the data to the output stringIO.
for path in paths:
    printData(path)

print "Written data to matrices.dat file."

outputFilename = "matrices.dat"
outputFile = open(outputFilename, 'w')
print >>outputFile, output.getvalue()[:-1]
outputFile.close()

print "Dumping panels from first pickle image..."
dumpPanels(paths)

# Now we dump image data for cppxfel.
threads = []
thread_count = int(os.getenv('NSLOTS', 4))
image_num = len(paths)
images_per_thread = image_num / thread_count

if (len(paths) == 0):
    print "No images with matrices provided."
    exit()

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

indexTxt = StringIO.StringIO()

print >> indexTxt, "ORIENTATION_MATRIX_LIST matrices.dat"
print >> indexTxt, "NEW_MATRIX_LIST indexed.dat"

print >> indexTxt, "# Be sure to set the UNIT_CELL and SPACE_GROUP for indexing. cppxfel cannot index without this knowledge."
print >> indexTxt, "SPACE_GROUP", spacegroup

if len(unit_cell_dimensions):
    print >> indexTxt, "UNIT_CELL",
    for dimension in unit_cell_dimensions:
        print >> indexTxt, dimension,
print >> indexTxt, "\n"

print >> indexTxt, "MM_PER_PIXEL", pixelSize[0]
print >> indexTxt, "BEAM_CENTRE", centre[0], centre[1]
print >> indexTxt, "DETECTOR_DISTANCE", distance
print >> indexTxt, "INTEGRATION_WAVELENGTH", wavelength, "\n"
print >> indexTxt, "PANEL_LIST panels.txt"
print >> indexTxt, "METROLOGY_SEARCH_SIZE 2\n"
print >> indexTxt, "# If your crystal is highly mosaic or the detector is quite far back you may need to increase the padding values."
print >> indexTxt, "SHOEBOX_FOREGROUND_PADDING 1"
print >> indexTxt, "SHOEBOX_NEITHER_PADDING 2"
print >> indexTxt, "SHOEBOX_BACKGROUND_PADDING 3\n"
print >> indexTxt, "# If you see too many spots, increase the intensity threshold.\nINTENSITY_THRESHOLD 12"
print >> indexTxt, "ABSOLUTE_INTENSITY OFF\n"
print >> indexTxt, "OVER_PRED_BANDWIDTH 0.07\n"
print >> indexTxt, "REFINE_ORIENTATIONS ON"
print >> indexTxt, "ROUGH_CALCULATION ON\n"

print >> indexTxt, "# Specifies maximum multiple lattices to index in total"
print >> indexTxt, "SOLUTION_ATTEMPTS 1\n"

print >> indexTxt, "# Maximum reciprocal distance from spot to spot to consider for analysis."
print >> indexTxt, "# A maximum reciprocal distance of 0.1 would be equivalent separation"
print >> indexTxt, "# between the beam centre and the 10 Angstrom resolution ring."
print >> indexTxt, "MAX_RECIPROCAL_DISTANCE 0.15\n"

print >> indexTxt, "# Initial rlp size: used to determine the tolerances for the vector lengths in the crystal."
print >> indexTxt, "# For a 1 micron crystal with no mosaicity, the initial rlp size is 0.0001 Ang^-1 (i.e.,"
print >> indexTxt, "# 1 / 10000 Ang). To be more strict for indexing, lower this number; to be less strict increase it."
print >> indexTxt, "INITIAL_RLP_SIZE 0.0001\n"

print >> indexTxt, "# If you wish to see more verbose output, change to 1 (moderate), or 2 (debug, usually too much)."
print >> indexTxt, "VERBOSITY_LEVEL 0\n"
print >> indexTxt, "COMMANDS\n"
print >> indexTxt, "INDEX"

outputFilename = "index.txt"

if not os.path.isfile(outputFilename):
    outputFile = open(outputFilename, 'w')
    print >>outputFile, indexTxt.getvalue()
    outputFile.close()
    print "New template input file index.txt"
    print "Please check your target unit cell and space group"


integrateTxt = StringIO.StringIO()

print >> integrateTxt, "ORIENTATION_MATRIX_LIST matrices.dat"
print >> integrateTxt, "# Enable this line to take results from cppxfel indexing (as opposed to DIALS):"
print >> integrateTxt, "# ORIENTATION_MATRIX_LIST integrate-indexed.dat"

print >> integrateTxt, "NEW_MATRIX_LIST orientations.dat\n"

print >> indexTxt, "# Be sure to set the SPACE_GROUP for integration."
print >> integrateTxt, "SPACE_GROUP", spacegroup

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
print >> integrateTxt, "# If your crystal is highly mosaic or the detector is quite far back you may need to increase the padding values."
print >> integrateTxt, "SHOEBOX_FOREGROUND_PADDING 1"
print >> integrateTxt, "SHOEBOX_NEITHER_PADDING 2"
print >> integrateTxt, "SHOEBOX_BACKGROUND_PADDING 3\n"
print >> integrateTxt, "# If you see too many spots, increase the intensity threshold."
print >> integrateTxt, "INTENSITY_THRESHOLD 12"
print >> integrateTxt, "ABSOLUTE_INTENSITY OFF\n"
print >> integrateTxt, "REFINE_ORIENTATIONS ON\n"
print >> integrateTxt, "# Enable the following line to be a rougher, quicker calculation"
print >> integrateTxt, "# ROUGH_CALCULATION ON\n"
print >> integrateTxt, "VERBOSITY_LEVEL 0\n"
print >> integrateTxt, "COMMANDS\n"
print >> integrateTxt, "INTEGRATE"

outputFilename = "integrate.txt"

if not os.path.isfile(outputFilename):
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

if not os.path.isfile(outputFilename):
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

if not os.path.isfile(outputFilename):
    outputFile = open(outputFilename, 'w')
    print >>outputFile, mergeTxt.getvalue()
    outputFile.close()
    print "New template input file merge.txt"

if distance == 0:
    print "You MUST change the experimental parameters in index.txt and integrate.txt"


print "Now what?"
print "********************"
print "***** INDEXING *****"
print "********************"
print
print "To index, change the space group and unit cell dimensions in index.txt."
print "After indexing, alter integrate.txt to accept the alternative value for"
print "ORIENTATION_MATRIX_LIST (uncomment the line)."
print
print "Command to index: cppxfel.run -i index.txt"
print
print "********************"
print "*** INTEGRATION ****"
print "********************"
print
print "To integrate orientations from indexing using DIALS, simply skip index.txt."
print "CHECK the unit cell and space group listed in integrate.txt, and then"
print "run the integrate.txt input file. This includes initial orientation matrix"
print "refinement."
print
print "Command to integrate: cppxfel.run -i integrate.txt"
print
print "... but first, wait for image dumps:"
