#!/usr/bin/python

import subprocess

class alignment_runner:
    
    # Defaults
    __SEQALTO_DEFAULT__ = "seqalto"
    __NUM_THREADS_DEFAULT__ = 2
    __ALIGNER_RESULTS_DEFAULT__ = "results.sam"

    def __init__(self):
        
        self.seqalto = self.__SEQALTO_DEFAULT__
        self.num_threads = self.__NUM_THREADS_DEFAULT__
        self.results_file = self.__ALIGNER_RESULTS_DEFAULT__
        
        print "Alignment runner initialized."

    def set_seqalto(self, seqalto_name):
        '''
            Sets the SeqAlto function caller, in case seqalto-basic or
            seqalto-fast is desired.
        '''
        
        self.seqalto = seqalto_name

    def set_num_threads(self, num_threads):
        '''
            Sets the number of threads to use.
        '''
        
        self.num_threads = num_threads
        
    def set_results_file(self, results_file):
        '''
            Sets the output alignment results filename.
        '''
        
        self.results_file = results_file

    def run_alignment(self, index_file, reads1_file, reads2_file):
        '''
            Runs SeqAlto on the provided reads files against the index.
            Returns the filename of the resulting .sam file.
        '''
    
        aligner_command = [self.seqalto, "align"]
        aligner_command.extend(["--idx", index_file])
        aligner_command.extend(["-1", reads1_file])
        aligner_command.extend(["-2", reads2_file])
        aligner_command.extend(["-p", str(self.num_threads)])

        print "Executing the command '%s'." % ' '.join(aligner_command)
        
        with open(self.results_file,"wb") as out:
            process = subprocess.Popen(aligner_command, stdout=out)
            
            try: count = (int(process.communicate()[0][:-1]))
            except: count = 0

        return self.results_file
    
    def get_num_reads(self, reads_file):
        '''
            Gets the approximate number of reads in the reads file.
        '''
        
        get_reads_command = ["wc", "-l", reads_file]
        
        p = subprocess.Popen(get_reads_command, stdout = subprocess.PIPE)
        total_lines = p.communicate()[0].split(" ")[0]
            
        try: count = (int(process.communicate()[0][:-1]))
        except: count = 0
        
        # Divide by 4 since each read has 4 lines associated with it
        # in a typical .fa file
        num_reads = int(total_lines) / 4
        
        return num_reads
        
        