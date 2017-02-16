#!/usr/bin/python

'''
File: bisig-prepare.py
Created: September 17, 2015
Authors: Paul Kowalski <paulkowa@buffalo.edu>
         Dhanasekar Karuppasamy <dhanasek@buffalo.edu>
Copyright (c) 2015-2016 Paul Kowalski, Dhanasekar Karuppasamy

Distributed under the MIT License.
See accompanying file LICENSE_MIT.txt.
This file is part of BiSiG.
'''

import os, sys
import re
import time, timeit
from Bio import SeqIO
from optparse import Option, OptionParser, OptionValueError
from os.path import dirname

inputMessage = "read input from this file/directory - must not contain file extensions"
outputMessage =  "write output to files with this prefix"
aminoMessage = "input is an amino acid sequence (default False)"
lengthMessage = "remove sequences shorter than this length (default 0)"
cleanMessage = "remove sequences with missing bases (default False)"
usageMessage = "--input name --output name [options...]" + "\n\nOptions:" + "\n  --input name" + "\n  -i name\t\t" + inputMessage + "\n  --output name" + "\n  -o name\t\t" + outputMessage + "\n  --amino {True|False}" + "\n  -a {True|False}\t" + aminoMessage + "\n  --length size" + "\n  -l size\t\t" + lengthMessage + "\n  --clean {True|False}" + "\n  -c {True|False}\t" + cleanMessage

def setUp():
    '''
    Handle command line arg options
    '''
    print "BiSiG Preparation Tool"
    print "Copyright (c) 2015-2016 SCoRe Group"
    print
    usagemessage = "Usage: %prog " + usageMessage
    parser = OptionParser(usage=usagemessage)

    parser.add_option("-i", "--input", \
        action="store", \
        type="string", \
        dest="input", \
        default=None, \
        help=inputMessage, \
        metavar="{FILE|DIR}")
    parser.add_option("-o", "--output", \
        action="store", \
        type="string", \
        dest="output", \
        default=None, \
        help=outputMessage, \
        metavar="FILE")
    parser.add_option("-a", "--amino", \
        action="store_false", \
        dest="dna", \
        default=True, \
        help=aminoMessage)
    parser.add_option("-l", "--length", \
        action="store", \
        type="int", \
        dest="length", \
        default=0, \
        help=lengthMessage)
    parser.add_option("-c", "--clean", \
        action="store_true", \
        dest="clean", \
        default=False, \
        help=cleanMessage)
    (options, sys.argv) = parser.parse_args()

    # Print error message if required options are not provided
    if (not options.input or not options.output or options.output.find('.') != -1 or options.output[-1:] == "/" or len(sys.argv) == 1):
        raise OptionValueError("")
    return options

def run(options):
    '''
    Run and print all information
    '''
    # Get list of files to be processed
    stats = Stats()
    fileList = []
    print "scanning " + options.input + " for input files..."
    files = checkDir(options, fileList)
    print "found " + str(files) + " input file(s)\nextracting sequences..."
    # Parse files
    start_time = timeit.default_timer()
    parse(options, fileList, stats)
    elapsed = timeit.default_timer() - start_time

    # Finish
    print "valid " + str(stats.valid) + " out of " + \
    str(stats.count) + " sequences\nwriting output files...\n" + \
    "shortest sequence: " + str(stats.short) + "\nlongest sequence: " + \
    str(stats.long) + "\naverage sequence: " + str(stats.avg) + \
    "\ntime: " + str("%.2f" % elapsed) + " seconds\ndone!"

def checkDir(options, filelist):
    '''
    Check input directory for all fasta / fsa files
    '''
    inputFiles = 0
    if options.input[-1:] == '/':
        for files in os.listdir(options.input):
            if files.endswith(".fa") or files.endswith(".fsa") or files.endswith(".fasta"):
                filelist.append(str(files))
                inputFiles += 1
    else:
        filelist.append(options.input)
        inputFiles += 1
    return inputFiles

def setOutput(options):
    if options.output.find('/') != -1:
        index = options.output.rfind('/')
        path = options.output[0:index + 1]
        output = options.output[index + 1:]
        os.chdir(path)
        return output
    return options.output


class Stats(object):
    '''
    Stores all statistical information for run
    '''
    def __init__(self):
        self.count = 0
        self.valid = 0
        self.avg = -1
        self.short = -1
        self.long = -1

    def update(self, length, valid):
        '''
        Update values in class
        '''
        self.shorter(length)
        self.longer(length)
        self.average(length, self.count)
        self.valid = valid
        self.count += 1

    def shorter(self, length):
        '''
        Check Shortest
        '''
        if self.short == -1:
            self.short = length
        elif length < self.short:
            self.short = length

    def longer(self, length):
        '''
        Check Longest
        '''
        if self.long == -1:
            self.long = length
        elif self.long < length:
            self.long = length

    def average(self, length, count):
        '''
        Compute rolling average
        '''
        if count == 0:
            self.avg = length
        else:
            self.avg = (self.avg + ((length - self.avg) / (count + 1)))

def parse(options, fileList, stats):
    '''
    Parse the input file and print to output file cleaned sequences
    '''
    valid = 0

    # Open output files
    output = setOutput(options)
    cleanFile = open(output + ".btxt", "w+")
    mapFile = open(output + ".bmap", "w+")
    removedSeqs = open(output + ".brm", "w+")

    # Iterate through input and write to output files
    for files in fileList:
        if options.input[-1:] == '/':
            inFile = open(options.input + files, 'rU')
        else:
            inFile = open(files, 'rU')

        for record in SeqIO.parse(inFile, "fasta"):
            # Process sequence and write to output files
            # Remove sequences < length
            # If clean remove all sequences containing invalid chars
            if options.clean:
                if len(record.seq) > options.length and \
                checkAlphabet(str(record.seq).upper(), options):
                    cleanFile.write(str(stats.valid) + '\t' + \
                        str(record.seq).upper() + '\n')
                    mapFile.write(str(stats.valid) + '\t' + \
                        str(record.description) + '\n')
                    valid += 1
                else:
                    removedSeqs.write(str(stats.count) + '\t' + \
                        str(record.description) + '\n')
            # If !clean replace all invalid chars with A
            else:
                if len(record.seq) > options.length:
                    if checkAlphabet(str(record.seq).upper(), options):
                        cleanFile.write(str(stats.valid) + '\t' + \
                            str(record.seq).upper() + '\n')
                        mapFile.write(str(stats.valid) + '\t' + \
                            str(record.description) + '\n')
                    else:
                        out = fixSeq(str(record.seq).upper(), options)
                        cleanFile.write(str(stats.valid) + '\t' + out + '\n')
                        mapFile.write(str(stats.valid) + '\t' + \
                            str(record.description) + '\n')
                    valid += 1
                else:
                    removedSeqs.write(str(stats.count) + '\t' + \
                        str(record.description) + '\n')
            # Record longest and shorest sequences
            stats.update(len(record.seq), valid)
        inFile.close()

    cleanFile.close()
    mapFile.close()
    removedSeqs.close()

def checkAlphabet(seq, options):
    '''
    Checks if sequence contains only valid characters
    '''
    # Check if valid DNA sequence
    if options.dna:
        if re.match("^[ACTG]*$", seq):
            return True

    # Check if valid protein sequence
    else:
        if re.match("^[ACDEFGHIKLMNPQRSTVWY]*$", seq):
            return True
    return False

def fixSeq(seq, options):
    '''
    Replace invalid characters in a sequence with A
    '''
    for letter in seq:
        if checkAlphabet(letter, options) == False:
            seq = seq.replace(letter, 'A')
    return seq

def main():
    '''
    Main
    '''
    try:
        run(setUp())
    except OptionValueError:
        print "Usage: bisig-prepare.py " + usageMessage
        print
        print '! error: incorrect command line arguments'
if __name__ == "__main__":
    main()
