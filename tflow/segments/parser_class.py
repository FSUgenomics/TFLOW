#TFLOW Component: OutputParser Base Class for Analyzing Job output
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

from time import sleep
import os.path
import sys
import subprocess
from ..util import print_except, print_exit

class OutputParser():
    def __init__(self):
        self.milestones = []
        self.out_file = ''
        self.job_type = 'Unspecified'
        self.terminal_flags = []
        self.failure_flags = []
        self.sleep_time = 5
        self.next_milestone_index = 0 
        self.current_milestone = 'Not Started'
        self.last_line = 'Tracking Not Started'
        self.last_line_index = 0
        self.shell_tracking = True
        self.tail_length = 15
        self.done_file_name = None
        self.running = True
        self.set_local_defaults()

    def set(self, settings={}, **kwargs):
        settings.update(kwargs)        
        for key in settings:
            if hasattr(self, key):
                setattr(self, key, settings[key])
            else:
                print_except('Unknown Key: %s Given to Parser' % key)
                

    def set_local_defaults(self):
        pass

    def track(self, loud=False):
        while self.running:
            if self.check_updated():
                self.running = self.check(loud)
            if self.running:
                sleep(self.sleep_time)

    def output_exists(self):
        return os.path.isfile(self.out_file)

    def read(self):
        if self.output_exists():
            out_file = open(self.out_file, 'r')
            for line in out_file:
                if '\r' in line:
                    print line.split('\r')[-1].rstrip()
                else:
                    print line.rstrip()
            out_file.close()
            
    def read_or_notify(self):
        if self.output_exists():
            self.read()
        else:
            print 'Output File: %s Does Not Exist.' % self.out_file

    #Add Last Milestone as Termination Flag
    def prepare_terminal_flags(self):
        if self.milestones:
            self.terminal_flags.append(self.milestones[-1])

    def check_terminal(self, line):
        for terminal_flag in self.terminal_flags:
            if terminal_flag in line:
                return True
        return False

    def check_failure(self, line):
        for failure_flag in self.failure_flags:
            if failure_flag in line:
                return True
        return False

    def completion_percent(self):
        return str(int((self.next_milestone_index+1)/float(len(self.milestones)) * 100)
                         ).zfill(2) + '% Completion '

    def check_completion(self, failure_exit=True):
        # If No Output File, Not Done.
        if not self.output_exists():
            return False

        #Get Last Output:
        if self.shell_tracking:
            output = subprocess.check_output(['tail', '-n', str(self.tail_length), self.out_file])
            output = output.splitlines()

        else:
            out_file = open(self.out_file, 'r')
            output = out_file.readlines()
            out_file.close()
        
            if len(output) > self.tail_length:
                output = output[(-1*self.tail_length):] 

        # If Empty Output File, Not Done.
        if not output:
            return False

        # If Failure Flag Found, Report and Terminate Search.
        for line in reversed(output):
            if self.check_failure(line) and failure_exit:
                print '%s Job has failed, Printing Final Output:' % self.job_type
                print ' --- '
                for line in output:
                    print line.rstrip()
                print ' --- '
                print '(To Retry, Delete Output File: %s )' % self.out_file
                print ''
                print_exit('', 0)

        # If Terminal Flag Found, Completed.
        self.prepare_terminal_flags()
        for line in reversed(output):
            if self.check_terminal(line):
                return True

        # No Terminal Flags Found, Not Completed.
        return False

    def check_updated(self):
        if not self.output_exists():
            print_except('Output File: %s Not Found!' % self.out_file)

        if self.shell_tracking:
            last_line = subprocess.check_output(['tail', '-n', '1', self.out_file])

        else:
            with open(self.out_file, 'r') as out_file:
                for line in out_file:
                    pass
                last_line= line
        return (last_line.strip() != self.last_line.strip())

    #Returns True if Not Yet Terminated
    def check(self, all_print=True):
        still_running = True
        if self.output_exists():
            out_file = open(self.out_file, 'r')
            output = out_file.readlines()
            out_file.close()
        else:
            print_except('Output File: %s Not Found!' % self.out_file)

        for line in output[self.last_line_index:]:
            if '\r' in line:
                print line.split('\r')[-1].rstrip()
            else:
                print line.rstrip()
 
            if self.check_terminal(line):
                still_running = False
                     
            if self.milestones[self.next_milestone_index] in line and still_running:
                self.milestone = self.milestones[self.next_milestone_index]
                self.next_milestone_index += 1
                
                if self.next_milestone_index == len(self.milestones):
                    still_running = False

            self.last_line = line.strip()
            self.last_line_index += 1
                
        return still_running

