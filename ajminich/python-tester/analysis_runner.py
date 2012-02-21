#!/usr/bin/python

import subprocess
import os
from datetime import datetime

class analysis_runner:

    # Defaults
    __SEQALTO_DEFAULT__ = "seqalto"
    __WIGGLE_DEFAULT__ = 21
    __MIN_MAPQ_DEFAULT__ = 1
    __MIN_LEN_DEFAULT__ = 1
    __ERROR_FILE_DEFAULT__ = "errorlog"
    __REPORT_DIRECTORY_DEFAULT__ = "results/"
    __RESULT_FILE_DEFAULT__ = "results.csv"
    __SEPARATOR__ = ','

    def __init__(self):
        
        self.seqalto = self.__SEQALTO_DEFAULT__
        self.wiggle = self.__WIGGLE_DEFAULT__
        self.min_mapq = self.__MIN_MAPQ_DEFAULT__
        self.min_len = self.__MIN_LEN_DEFAULT__
        self.error_file = self.__ERROR_FILE_DEFAULT__
        self.report_directory = self.__REPORT_DIRECTORY_DEFAULT__
        self.result_file = self.report_directory + '/' + self.__RESULT_FILE_DEFAULT__
        
        print "Analysis runner initialized."
        
    def run_analyzer(self, map_files, map_names, aligner_results_file, num_reads):
        '''
            Runs SeqAlto pair-checking analysis on the provided alignment
            against the original simulated map file.
        '''
        
        # Create a directory for the reports
        time_string = datetime.now().strftime("%Y-%m-%d %H:%M")
        subfolder = "Results - " + time_string
        report_folder = self.report_directory + "/" + subfolder
        
        if not os.path.exists(report_folder):
            os.makedirs(report_folder)
        
        reports = []
        num_map_files = len(map_files)
        
        # Generate a report for each map file
        for index in range(num_map_files):
            
            map_file = map_files[index]
            map_name = map_names[index]
            
            num_reads_per_map_file = num_reads / num_map_files
            
            report = self.get_report(map_file, map_name, aligner_results_file, num_reads_per_map_file)
        
            # Save the report to a file
            current_report_file = report_folder + "/" + map_name + "_results.txt"
            f = open(current_report_file, 'w')
            f.write(report)
            f.close()
            
            # Append the report to the current list
            reports.append(report)
            
        # Parse the reports
        parsed_reports = [self.parse_report(report) for report in reports]
        
        # Print the reports to a final result file
        f = open(self.result_file, 'w')
        f.write("Results: " + time_string)
        f.write('\n\n' + self.__SEPARATOR__)
        
        for map_name in map_names:
            f.write(map_name + self.__SEPARATOR__)
        f.write('\n')
        
        for result in parsed_reports[0].keys():
            f.write(result)
            
            for parsed_report in parsed_reports:
                f.write(self.__SEPARATOR__ + parsed_report[result])
            
            f.write('\n')
        
    def get_report(self, map_file, map_name, aligner_results_file, num_reads):
        '''
            Gets the read-by-read comparison report for a single map file.
        '''
        
        map_results = map_name + ".sam"
            
        # cat results.sam | grep pa > results_pa.sam
        p1 = subprocess.Popen(["cat", aligner_results_file], stdout = subprocess.PIPE)
        
        with open(map_results, "wb") as out:
            p2 = subprocess.Popen(["grep", map_name], stdin = p1.stdout, stdout = out)
    
        p1.stdout.close()
        result = p2.communicate()[0]
        p1.wait()
    
        # seqalto check_pair randCEUma.fa.map results.sam 21 7200000 1 1 Kakashi > ma_results
        analyzer_command = [self.seqalto, "check_pair"]
        analyzer_command.append(map_file)
        analyzer_command.append(map_results)
        analyzer_command.append(str(self.wiggle))
        analyzer_command.append(str(num_reads))
        analyzer_command.append(str(self.min_mapq))
        analyzer_command.append(str(self.min_len))
        analyzer_command.append(self.error_file)

        print "Executing the command '%s'." % ' '.join(analyzer_command)
        
        p = subprocess.Popen(analyzer_command, stdout = subprocess.PIPE)
        report = p.communicate()[0]
        
        return report
    
    def parse_report(self, report_string):
        '''
            Returns a dictionary containing the results of a read-by-read
            alignment comparison. 
        '''
        
        results = [line.split(':') for line in report_string.split('\n')]
        
        report_info = {}
        for result in results:
            if len(result) >= 2:
                report_info[result[0].strip()] = result[1].strip()
            
        return report_info
        