#!/usr/bin/env libtbx.python

from multiprocessing import Process
import sys
from os.path import basename
import os
import scitbx_array_family_flex_ext
from subprocess import call
import StringIO
from cppxfel.command_line_parser import CommandLineParser

parser = CommandLineParser()
parser.setup()

import_phil = "import.options"
find_spots_phil = "find_spots.options"
index_phil = "index.options"

skip = 0
max = 0
index = True

try:
        shouldIndex = parser.valueForKey('index')
        if 'n' in shouldIndex:
                index = False
                print "Will not index"
except:
        pass

try:
        skip = int(parser.valueForKey('skip'))
except:
        pass

try:
        max = int(parser.valueForKey('max'))
except:
        pass

if len(sys.argv) == 1:
        print "Script to run DIALS on individual still shots (XFEL).\n"

        print "This script will run dials.import, dials.find_spots and dials.index"
        print "on image files provided individually. Any user-defined parameters"
        print "should be defined using the import.options, find_spots.options and"
        print "index.options files which will be tagged onto the end of the DIALS"
        print "invocation. Do not overwrite the input or output commands.\n"

        print "Supply the images on the command line using a Unix pattern to match"
        print "the appropriate files in the current directory. For example:\n"
        print "\tcppxfel.run_dials shot*.pickle\n"
        print "The number of cores will be fetched from the environment variable"
        print "$NSLOTS. If this is not set the default will be 4 cores."
        exit(1)

"""
Imports, spot-finds and indexes images using DIALS.

Args: images: a list of paths to image files in the current directory.

Returns: None

"""
def indexImages(images):
        global start, end
        jsons = []

        import_options = ""
        if os.path.isfile(import_phil):
                import_options_file = open(import_phil, 'r')
                import_options = import_options_file.read()

        print "Import options: ", import_options

        find_spots_options = ""
        if os.path.isfile(find_spots_phil):
                find_spots_options_file = open(find_spots_phil, 'r')
                find_spots_options = find_spots_options_file.read()

        find_spots_options = find_spots_options.replace("\n", "")
        print "Find spots options: ", find_spots_options

        index_options = ""
        if os.path.isfile(index_phil):
                index_options_file = open(index_phil, 'r')
                index_options = index_options_file.read()

        index_options = index_options.replace("\n", "")

        print "Index options: ", index_options

        for i in range(0, len(images)):
                filename = images[i]
                base = basename(filename)
                rootname = os.path.splitext(base)[0]

                experimentJson = "_" + rootname + "_experiments.json"
                if os.path.isfile(experimentJson):
                        print "Skipping", rootname, " - already indexed"
                        continue

                jsons.append(rootname)

                skip = os.path.isfile(rootname + ".json")

                if skip == False:
                        command = "dials.import " + filename + " output.datablock=" + rootname + ".json "
                        command += import_options
                        print "Executing command: ", command
                        os.system(command)
                else:
                        print "Skipping import for " + filename + ", already present"

                skip = os.path.isfile("_" + rootname + "_strong.pickle")

                if skip == False:
                        command = "dials.find_spots " + rootname + ".json output.reflections=_" + rootname + "_strong.pickle "
                        command += find_spots_options + " | tee " + rootname + ".find_spots.log"
                        print "Executing command:", command
                        os.system(command)
                else:
                        print "Skip spotfinding for " + rootname + ", already present"

                if not index:
                        continue

                command = "dials.index " + rootname + ".json _" + rootname + "_strong.pickle "
                command += " output.reflections=_" + rootname + "_indexed.pickle "
                command += "output.experiments=_" + rootname + "_experiments.json "
                command += index_options + " | tee " + rootname + ".index.log"
                print "Executing command:", command

                os.system(command)

"""
Function finds the filename stem from an experiment JSON file
excluding the extension.

Args: filename: string filename of JSON file.

Returns: filename stem.

"""
def rootnameStem(filename):
        return filename.replace("_experiments.json", "")

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
        print >>output, "image " + rootname_stem

        for experiment in experiments:
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
        rootname_stem = rootnameStem(filename)

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

        printExperiments(experiments, filename, output)

output = StringIO.StringIO()

threads = []

start = 0
end = len(parser.getLoose())

if (skip > 0):
        start = skip

if (max > 0):
        end = start + max

if start > len(parser.getLoose()):
        print "Start position is beyond end of supplied images."
        exit()

if end > len(parser.getLoose()):
        end = len(parser.getLoose())

myImages = parser.getLoose()[start:end]

thread_count = int(os.getenv('NSLOTS', 4))
image_num = len(myImages)
print "Total images:", image_num

images_per_thread = image_num / thread_count

for i in range(0, thread_count):
                min = int(i * images_per_thread)
                max = int((i + 1) * images_per_thread)

                print max - min, "images on this thread."

                thread = Process(target=indexImages, args=(myImages[min:max], ))
                threads.append(thread)
                thread.start()
