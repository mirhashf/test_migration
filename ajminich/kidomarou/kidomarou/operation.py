'''
Defines the types of operations that can be performed by the 
genome processing pipeline. Each edge in the Kidomarou graph 
is an operation.
'''

from .file_groups import *

class Operation:
    '''Defines a generic operation.
    
    Derived classes must implement a run() method.
    '''
    
    def __init__(self, name, input_file_group, output_file_group):
        
        self._name = name
        self._input_file_group = input_file_group
        self._output_file_group = output_file_group
        
    def get_name(self):
        return self._name
        
    def get_input_file_group(self):
        return self._input_file_group
    
    def get_output_file_group(self):
        return self._output_file_group

class BinaAlignmentOperation(Operation):
    
    NAME = "Read Alignment with Bina Aligner"
    
    def __init__(self):
        
        input_group = FastqGroup()
        output_group = AlignedReadsGroup()
        
        Operation.__init__(self, self.NAME, input_group, output_group)
    
    def run(self):
        '''Should run the read alignment.
        '''
        
        print "Not yet implemented!"