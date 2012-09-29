
'''
Enums to define:
- FileDivision (None, Coordinate, Lane)
- FileType (FASTA, SAM, BAM, VCF, GFF)
- SortOrder (None, Coordinate, Queryname)
'''



class FileGroup:
    '''Defines a file group. Each node in the Kidomarou  
    web is a FileGroup class, representing a set of data 
    or results files that can be generated from a series 
    of processes (edges).
    
    This object can be extended into derived classes as 
    necessary. However, the nodes of the graph must all 
    fundamentally be FileGroup objects.
    
    All classes based on this function must implement the 
    get_name() function to allow them to be stored in the 
    graph.
    '''
    
    _files = []
    
    def __init__(self, file_type, is_multiple_files):
        self._file_type = file_type
        self._is_multiple_files = is_multiple_files
    
    def get_file_type(self):
        return self._file_type
    
    def get_is_multiple_files(self):
        return self._is_multiple_files
    
    def get_files(self):
        return self._files
    
    def set_files(self, files):
        self._files = files
        
class FastqGroup(FileGroup):
    '''Represents a group of FASTQ reads files.'''
    
    FILE_TYPE = "fastq"
    
    def __init__(self, is_multiple_files):
        FileGroup.__init__(self, self.FILE_TYPE, is_multiple_files)
        
    def get_name(self):
        multi_name = ("multi" if self.get_is_multiple_files() else "single") 
        return self.FILE_TYPE + "_" + multi_name
    
class AlignedReadsGroup(FileGroup):
    '''Represents a group of aligned reads.'''
    
    def __init__(self, is_bam, is_multiple_files, division, sort_order):
        self._file_type = ("bam" if is_bam else "sam")
        self._division = division
        self._sort_order = sort_order
        
        FileGroup.__init__(self, self._file_type, is_multiple_files)
        
    def get_name(self):
        return (self.get_file_type() + "_" + self.get_division() + 
                "_" + self._get_sort_order())
    
    def get_division(self):
        return self._division
    
    def get_sort_order(self):
        return self._sort_order
        
