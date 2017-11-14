def printSpots(filenames):
    for filename in filenames:
        fin = open(filename, 'rb')
        strongPickle = pickle.load(fin)
        output = StringIO.StringIO()

        for spot in strongPickle:
            xyzobs = spot['xyzobs.px.value']
            x = xyzobs[0]
            y = xyzobs[1]
            print >>output, str(x) + "\t" + str(y)

        new_name = filename.replace('pickle', 'list')
        fout = open(new_name, 'w')
        fout.write(output.getvalue())
        fout.close()
        fin.close()

if __name__ == '__main__':
    import pickle
    import cppxfel
    import StringIO
    import sys

    majorOutput = StringIO.StringIO()
    printSpots(sys.argv[1:])
