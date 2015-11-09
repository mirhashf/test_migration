# adding 2 roman numerals without changing to arabian numbers

romanNumeralOrdering = "IVXLCDM"
subtractionMapper = {"IV":"IIII", "IX":"VIIII", "XL":"XXXX", "XC":"LXXXX", "CD":"CCCC", "CM":"DCCCC"}
additionMapper = {"IIIII":"V", "VV":"X", "XXXXX":"L", "LL":"C", "CCCCC":"D", "DD":"M"}


def addValues(first, second):

	new_first = removeSubtraction(first)
	new_second = removeSubtraction(second)

	# sort the numerals by their values
	sortedSum = ''.join(sorted(new_first+new_second, key=romanNumeralOrdering.index))[::-1]
	# add the numerals according to their mappings
	aggregatedSum = aggregateValue(sortedSum)
	# change the current values to subtractions, if they exist
	result = performSubtraction(aggregatedSum)
	
	print "the sum of %s and %s is : " %(first,second) + result

	return result
	


def aggregateValue(value):
	# replaace a group of smaller values with a larger one
	changed = True
	# have multiple passes to account for reductions that did not exist in first pass
	# i.e. "VIIIII" -> "VV" -> "X"
	while changed:
		changed = False
		for key, val in additionMapper.iteritems():
			if key in value:
				value = value.replace(key, val)
				changed = True
	return value

def removeSubtraction(value):
	# remove any subraction that exists in the roman numeral value
	# in order to add the two values together
	for key, val in subtractionMapper.iteritems():
		if key in value:
			value = value.replace(key, val)
	return value

def performSubtraction(value):
	# add in any subraction that exists in the roman numeral value
	for key, val in subtractionMapper.iteritems():
		if val in value:
			value = value.replace(val, key)
	return value

#  the function addValues returns the new roman numeral, however,
#  just for clarity there is a print statement to show the result
#  of the addition (could also store the result of a call to addValues)
addValues("CLI", "LIX")
addValues("CCCLXIX","DCCCXLV")
addValues("III","IX")
# Start Emre
addValues("I","II")
addValues("IV","VI")




