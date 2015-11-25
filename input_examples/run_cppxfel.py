#!/usr/bin/env libtbx.python

import cppxfel
parser = cppxfel.cppParser('integrate.txt')
parser.integrate()
