'''

    Method to add any two roman numerals together. 
    Highest letter that can be input is M (1000).
    Algorithm idea borrowed from: 
    http://turner.faculty.swau.edu/mathematics/materialslibrary/roman/. 

    Assumes input is valid uppercase roman numeral. 

'''

def roman_numeral_addition(first, second):
    first, second = replace_subtractives(first), replace_subtractives(second)

    combined = list(first + second)
    num_dict = combine_groups(combined)
    return dict_to_str(num_dict)

def replace_subtractives(numeral):
    return numeral.replace('CM', 'CCCCC').replace('CD', 'CCCC').replace('XC', 'LXXXX').replace('XL', 'XXXX').replace('IX', 'VIIII').replace('IV', 'IIII')

def sub_subtractives(numeral):
    return numeral.replace('CCCCC', 'CM').replace('CCCC', 'CD').replace('LXXXX', 'XC').replace('XXXX', 'XL').replace('VIIII', 'IX').replace('IIII', 'IV')

def combine_groups(vals):
    dict_vals = {'M': 0, 'D': 0, 'C': 0, 'L': 0, 'X': 0, 'V': 0, 'I': 0}
    for letter in vals:
        dict_vals[letter] += 1

    while dict_vals['I'] >= 5:
        dict_vals['I'] -= 5
        dict_vals['V'] += 1
    while dict_vals['V'] >= 2:
        dict_vals['V'] -= 2
        dict_vals['X'] += 1
    while dict_vals['X'] >= 5:
        dict_vals['X'] -= 5
        dict_vals['L'] += 1
    while dict_vals['L'] >= 2:
        dict_vals['L'] -= 2
        dict_vals['C'] += 1
    while dict_vals['C'] >= 5:
        dict_vals['C'] -= 5
        dict_vals['D'] += 1
    while dict_vals['D'] >= 2:
        dict_vals['D'] -= 2
        dict_vals['M'] += 1

    return dict_vals

def dict_to_str(num_dict):
    result = ''
    result += 'M' * num_dict['M']
    result += 'D' * num_dict['D']
    result += 'C' * num_dict['C']
    result += 'L' * num_dict['L']
    result += 'X' * num_dict['X']
    result += 'V' * num_dict['V']
    result += 'I' * num_dict['I']
    return sub_subtractives(result)

def mergesort(list1, list2):
    result = []

    while len(list1) > 0 and len(list2) > 0:
        if compare_letters(list1[0], list2[0]) < 0:
            result.append(list1.pop(0))
        else:
            result.append(list2.pop(0))
    while len(list1) > 0:
        result.append(list1.pop(0))
    while len(list2) > 0:
        result.append(list2.pop(0))
    return result

def compare_letters(num1, num2):
    if num1 == num2:
        return 0
    elif num1 == 'M':
        return -1
    elif num2 == 'M':
        return 1
    elif num1 == 'I':
        return 1
    elif num2 == 'I':
        return -1
    elif num1 == 'V':
        if num2 == 'I':
            return -1
        else:
            return 1
    elif num1 == 'X':
        if num2 == 'I' or num2 == 'V':
            return -1
        else:
            return 1
    elif num1 == 'L':
        if num2 == 'I' or num2 == 'V' or num2 == 'X':
            return -1
        else:
            return 1 
    elif num1 == 'C':
        if num2 == 'I' or num2 == 'V' or num2 == 'X' or num2 == 'L':
            return -1
        else:
            return 1 
    elif num1 == 'D':
        if num2 == 'I' or num2 == 'V' or num2 == 'X' or num2 == 'C':
            return -1
        else:
            return 1 

def confirm_valid_arg(numeral):
    values = set(list(numeral))
    for val in values:
        if val not in (['M', 'D', 'C', 'L', 'X', 'V', 'I']):
            print 'Not a valid argument.'
            print 'Roman numerals are composed of the following characters: M, D, C, L, X, V, I'
            return False
        
    return True

print 'EXAMPLE: Addition of %s and %s is: %s' % ('CCCLXIX', 'DCCCXLV', roman_numeral_addition('CCCLXIX', 'DCCCXLV'))

if __name__ == '__main__':
    while True:
        print 'Please enter first numeral.'
        val1 = raw_input()
        if not confirm_valid_arg(val1):
            continue
        print 'Please enter second numeral.'
        val2 = raw_input()
        if not confirm_valid_arg(val2):
            continue
        print 'Addition of %s and %s is: %s' % (val1, val2, roman_numeral_addition(val1, val2))

