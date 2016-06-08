# Quick command line parser

import sys

# takes string "key"
class CommandLineParser:
	dict = {}
	loose = []
	
	def setup(self):
		for arg in sys.argv[1:]:
			key_value = arg.split("=")
		
			if (len(key_value) < 2):
				self.loose.append(key_value[0])
				continue
			else:
				values = key_value[1].split(",")
				
				if (len(values) == 1):
					self.dict[key_value[0]] = key_value[1]
				else:
					self.dict[key_value[0]] = values
	
	def valueForKey(self, key):
		if key in self.dict:
			return self.dict[key]
		else:
			return ""
			
	def empty(self):
		return (len(self.dict) == 0)
	
	def getLoose(self):
		return self.loose

