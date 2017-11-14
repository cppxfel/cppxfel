import pickle
import StringIO

def printSpots(filenames):
    for filename in filenames:
        fin = open(filename, 'rb')
        strongPickle = pickle.load(fin)
        output = StringIO.StringIO()

        for spot in strongPickle['xyzobs.px.value']:
            x = spot[0]
            y = spot[1]
            print >>output, str(x) + "\t" + str(y)

        new_name = filename.replace('pickle', 'list')
        fout = open(new_name, 'w')
        fout.write(output.getvalue())
        fout.close()
        fin.close()

# if __name__ == '__main__':
#       import pickle
#       import cppxfel
#
#
#
#       majorOutput = StringIO.StringIO()
#       printSpots(sys.argv[1:])
