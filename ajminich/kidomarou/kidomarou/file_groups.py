'''
Defines the types of file groups that can exist in the genome
processing pipeline. Each node in the Kidomarou graph is a 
file group.
'''

class FileGroup:
    '''Defines a file group. Each node in the Kidomarou  
    graph is a FileGroup class, representing a set of data 
    or results files that can be generated from a series 
    of processes (edges).
    
    This object can be extended into derived classes as 
    necessary. However, the nodes of the graph must all 
    fundamentally be FileGroup objects.
    
    All classes based on this function must implement the 
    get_id() function to allow them to be stored in the 
    graph.
    '''
    
    # Types
    FASTQ = "fastq"
    BAM = "bam"
    SAM = "sam"
    VCF = "vcf"
    GFF = "gff"
    
    _files = []
    
    def __init__(self, name, file_type, is_multiple_files):
        self._name = name
        self._file_type = file_type
        self._is_multiple_files = is_multiple_files
    
    def get_name(self):
        return self._name
    
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
    
    def __init__(self, name, is_multiple_files):
        FileGroup.__init__(self, name, FileGroup.FASTQ, is_multiple_files)
        
    def get_id(self):
        multi_name = ("multi" if self.get_is_multiple_files() else "single") 
        return FileGroup.FASTQ + "_" + multi_name
    
class AlignedReadsGroup(FileGroup):
    '''Represents a group of aligned reads.'''

    _is_realigned = False
    _is_recalibrated = False
    
    def __init__(self, name, is_bam, is_multiple_files, division, sort_order):
        self._file_type = (FileGroup.BAM if is_bam else FileGroup.SAM)
        self._division = division
        self._sort_order = sort_order
        
        FileGroup.__init__(self, name, self._file_type, is_multiple_files)
        
    def get_id(self):
        name_base = (self.get_file_type() + "_" + self.get_division() + 
                "_" + self.get_sort_order()) 
        
        realigned_state = ("_realigned" if self._is_realigned else "")
        recal_state = ("_recal" if self._is_recalibrated else "")
        
        return name_base + realigned_state + recal_state
    
    def get_division(self):
        return self._division
    
    def get_sort_order(self):
        return self._sort_order
        
    def set_is_realigned(self, is_realigned):
        self._is_realigned = is_realigned
        
    def get_is_realigned(self):
        return self._is_realigned
    
    def set_is_recalibrated(self, is_recalibrated):
        self._is_recalibrated = is_recalibrated
        
    def get_is_recalibrated(self):
        return self._is_recalibrated
    
class VariantCallsGroup(FileGroup):
    '''Represents a group of called variants.'''

    _is_indel_left_aligned = True
    _is_recalibrated = False
    
    def __init__(self, name, is_vcf, is_multiple_files, division):
        self._file_type = (FileGroup.VCF if is_vcf else FileGroup.GFF)
        self._division = division
        
        FileGroup.__init__(self, name, self._file_type, is_multiple_files)
        
    def get_id(self):
        name_base = self.get_file_type() + "_" + self.get_division() 
        
        realigned_state = ("_left_aligned" if self._is_indel_left_aligned else "")
        recal_state = ("_recal" if self._is_recalibrated else "")
        
        return name_base + realigned_state + recal_state
    
    def get_division(self):
        return self._division
    
    def get_sort_order(self):
        return self._sort_order
        
    def set_is_indel_left_aligned(self, is_indel_left_aligned):
        self._is_indel_left_aligned = is_indel_left_aligned
        
    def get_is_indel_left_aligned(self):
        return self._is_indel_left_aligned
    
    def set_is_recalibrated(self, is_recalibrated):
        self._is_recalibrated = is_recalibrated
        
    def get_is_recalibrated(self):
        return self._is_recalibrated
