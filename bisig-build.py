#!/usr/bin/python

from __future__ import print_function
import sys
import random
import struct
import hashlib
import os, sys
import itertools
import numpy as np
import operator as op
from ctypes import c_uint64
from os.path import dirname, isdir
from optparse import Option, OptionParser
from pyspark import SparkContext

#Begin GUI Section
kmerlengthMessage = "use kmers of this size (default 16)"
moduloMessage = "use this mod value in sketching (default 20)"
iterationsMessage = "limit the number of sketching iterations to this (default 5)"
cmaxMessage = "use this limit to mark frequent kmers (default 15000)"
jminMessage = "use this limit to extract candidate pairs (default 85)"
partitionsMessage = "run on this many partitions (default 1024)"
tmpfileMessage = "use this temporary file during program execution (default 'tmp.py')"
inputMessage = "read input from this file (default '/home/in.btxt')"
outputMessage =  "write output to files with this prefix (default '/home/out')"
smhasherMessage = "the location of the smhasher library (default './smhasher-install')"
usageMessage = "--input name --output name [options...]" \
+ "\n\nOptions:" \
+ "\n  --input name" + "\n  -i name\t\t" + inputMessage \
+ "\n  --output name" + "\n  -o name\t\t" + outputMessage \
+ "\n  --smhasher name" + "\n  -s name\t\t" + smhasherMessage \
+ "\n  --temp name" + "\n  -t name\t\t" + tmpfileMessage \
+ "\n  --kmer size" + "\n  -k size\t\t" + kmerlengthMessage \
+ "\n  --modulo number" + "\n  -m number\t\t" + moduloMessage \
+ "\n  --iter number" + "\n  -e number\t\t" + iterationsMessage \
+ "\n  --cmax number" + "\n  -c number\t\t" + cmaxMessage \
+ "\n  --jmin number" + "\n  -j number\t\t" + jminMessage \
+ "\n  --partition number" + "\n  -p number\t\t" + partitionsMessage

cmdLineError = False
def errorMessage(errorCode):
    '''
    Generates and prints an error message based on the specified error.
    '''
    cmdLineError = True
    print("Recommended: redirect standard error")
    print("Usage: spark-submit bisig-build.py " + usageMessage)
    print()
    if(errorCode == 1):#sc = SparkContext()
        print('! error: SparkContext could not be gathered')
        print('! ensure that you are running through spark-submit and that start-master has been run')
    elif(errorCode == 2):#Input
        print ('! error: input filepath does not exist')
    elif(errorCode == 3):#Output
        print ('! error: output filepath is invalid or already exists')
    elif(errorCode == 4):#Kmer Length
        print ('! error: kmer length is too short')
    elif(errorCode == 5):#Modulo
        print ('! error: modulo value is too low')
    elif(errorCode == 6):#Iterations
        print ('! error: number of iterations is invalid')
    elif(errorCode == 7):#Cmax
        print ('! error: cmax is too small')
    elif(errorCode == 8):#Jmin
        print ('! error: jmin is too small')
    elif(errorCode == 9):#Partitions
        print ('! error: number of partitions is invalid')
    elif(errorCode == 10):#Smhasher Location
        print ('! error: smhasher location is invalid')
    else:
        print ('! error: incorrect command line arguments')

#Required for Spark
try:
    sc = SparkContext()
except:
    errorMessage(1)

#End of GUI code at start of file

class options:
    kmer_length = 16
    modulo = 20
    iterations = 5
    cmax = 15000
    jmin = 85
    partitions = 1024
    tmp_file = "tmp.py"
    ip = "/home/in.btxt"
    output = "/home/out"
    smhasher_loc = "./smhasher-install"

def gen_kmers(input, k):
    '''
    Iterator to generate KMERS for a given string.
    '''
    for i in xrange(len(input) - k + 1):
        yield (str(input[i:k + i]), i)

def map_sketch(seq, kmer_length, modulo, i):
    '''
    Generates <hash, <seq_id, kmer_position, kmer_count>> for every Kmer in the sequence.
    Uses a smhasher python binding to generate the hashes to be in 1:1 mapping with Elastic.
    '''
    sys.path.append(dirname(options.smhasher_loc + "/lib/python/"))
    from smhasher import murmur2
    seq_id, seq_txt = seq.split('\t')
    kmer_list = []
    for kmer_str, kmer_pos in gen_kmers(seq_txt, kmer_length):
        kmer_hash = murmur2(kmer_str,2147483647)
        if(kmer_hash%modulo==i):
            kmer_list.append((kmer_hash, kmer_pos))
        else:
            continue
    kmer_count = len(kmer_list)
    for kmer_hash, kmer_pos in kmer_list:
        yield (kmer_hash, [(int(seq_id), kmer_pos, kmer_count)])

def combine_pairs(iterator):
    '''
    Old one. Not used anymore. Look @ combine_pairs_fast which considers rectangles & triangles.
    Uses itertools combinations (https://docs.python.org/2.7/library/itertools.html#itertools.combinations) 
    to create sketch pairs. 
    Emits <<seqid1, seqid2, <min length of seq1, seq2>, <number of matching seq>>
    '''
    for x in iterator:
        for i in itertools.combinations(x[1], 2):
            sketchi, sketchj = i
            if(sketchi[0] < sketchj[0]):
                si = sketchi
                sj = sketchj
            elif(sketchj[0] < sketchi[0]):
               si = sketchj
               sj = sketchi
            yield ((si[0], sj[0], min(si[2],sj[2])), (si[1], sj[1]))

def build_aux(sketch):
    '''
    Builds auxillary list <seqid, [hash]>
    '''
    for seqs in sketch[1]:
        yield (seqs[0], [sketch[0]])

def update_counts(x):
    '''
    Updates counts on the auxillary list. This was the last piece of code written for 1:1 mapping . Paul knows this logic.
    '''
    for pair in x:
        if not (auxList.value.get(pair[0][0]) is None) and not (auxList.value.get(pair[0][1]) is None):
            sketchi = auxList.value.get(pair[0][0])
            sketchj = auxList.value.get(pair[0][1])
            yield ((pair[0][0], pair[0][1]), int(((float(pair[1]) + float(len(set(sketchi) & set(sketchj)))) / float(pair[0][2])) * 100.0))
        else:
            yield ((pair[0][0], pair[0][1]), int(((float(pair[1]) / float(pair[0][2])) * 100.0)))

def stats(iterator):
    '''
    Utility function to get stats about number of items in a partition.
    Use rdd.mapPartition(stats) to get the details of stats. 
    '''
    count = 0
    for i in iterator:
        count = count + 1
        pass
    yield count

def rand_part(i):
    '''
    Random number generator for partition allocation used for normalised distibution.
    '''
    return random.randint(0,options.partitions)

def split_list(lis, split=20):
    '''
    Splits a list in to nested sublist with 20 items each.
    '''
    lis_len = (len(lis)/split)+1 if (len(lis) % split != 0) else (len(lis)/split)
    return [lis[i*split:(i*split)+split] for i in range(lis_len)]

def break_long_ones(iterator):
    '''
    Breaks long list in to triangles and rectangles.
    1. Split the lists to nested list of 20 each.
    2. For each item in the nested list
        - emit ('s', list) -> Rectangles.
        - emit ('p', (list1, list2)) -> Traingles.
    
    https://docs.python.org/2.7/library/itertools.html#itertools.combinations takes care of pair wise iteration.
    '''
    for i in iterator:
        single = split_list(list(i[1]))
        for l in single:
            yield ('s', l)
        for l in itertools.combinations(single, 2):
            yield ('p', l)

def combine_pairs_fast(iterator):
    '''
    Combine pairs based on if they are triangles or rectangles.
    As a first step look @ the combine_pairs code above to understand how this has evolved.
    '''
    for x in iterator:
        if(x[0] == 's'):
            for sketchi, sketchj in itertools.combinations(x[1], 2):
                if(sketchi[0] < sketchj[0]):
                    si = sketchi
                    sj = sketchj
                elif(sketchj[0] <= sketchi[0]):
                   si = sketchj
                   sj = sketchi
                yield ((si[0], sj[0], min(si[2],sj[2])), (si[1], sj[1]))
        elif(x[0] == 'p'):
            for sketchi, sketchj in itertools.product(x[1][0], x[1][1]):
                if(sketchi[0] < sketchj[0]):
                    si = sketchi
                    sj = sketchj
                elif(sketchj[0] <= sketchi[0]):
                   si = sketchj
                   sj = sketchi
                yield ((si[0], sj[0], min(si[2],sj[2])), (si[1], sj[1]))

def iteration(options):
    '''
    Generates the meta code in tmp.py used for iterations and runs it using eval.
    A file called tmp.py is generated when this code is executed.
    
    Code generated for 1 iteration.
    
    # Load file
    btxtRDD = sc.textFile(options.ip, 1024)
    global broadcastAggAuxList, aggAuxList, candidateRDD, outRDD, auxList
    auxList = None
    
    
    # For iteration 0 Map Every KMER to <hash, <seq_id, kmer_position, kmer_count>>
    # Reduce it by hash to make it <Hash : <list>>
    sketchRDD0 = btxtRDD.flatMap(lambda s: (map_sketch(s, options.kmer_length, options.modulo,0))).reduceByKey(lambda a, b: (a + b)).cache()
    
    # Build Aux RDD by filtering ones which have long lists. Then build a list of <Seqid, Hash> then reduce them to make is <Seqid, <List of Hashes>> 
    auxRDD0 = sketchRDD0.filter(lambda v: len(v[1]) >= options.cmax).flatMap(lambda s: build_aux(s)).reduceByKey(lambda a, b: a + b)
    auxList = sc.broadcast(auxRDD0.collectAsMap())
    
    
    # Combine pairs.
    combinePairsRDD0 = sketchRDD0.filter(lambda v: len(v[1]) > 1 and len(v[1]) < options.cmax).mapPartitions(break_long_ones).partitionBy(options.partitions, partitionFunc=rand_part).mapPartitions(combine_pairs_fast)
    candidateRDD0 = combinePairsRDD0.aggregateByKey(0, lambda a, b: a + 1, lambda x, y: x + y).mapPartitions(update_counts).filter(lambda a: a[1] >= options.jmin)
    globals()['candidateRDD0'] = candidateRDD0
    
    
    # Union. Answer.
    candidateRDD=candidateRDD0.union(candidateRDD1).union(candidateRDD2)
    globals()['candidateRDD'] = candidateRDD
    candidateRDD.cache()
    
    TODO : Meta code is too cryptic. Find a better way for iteration.
    
    '''
    code = "inputSeqRDD = sc.textFile(options.ip, {partitions})".format(partitions = options.partitions)
    code = code + '\n'+ "global auxList" \
        + '\n' + "auxList = None"
    for i in range(options.iterations):
        code = code + '\n' + "sketchRDD{iter_number} = inputSeqRDD.flatMap(lambda s: (map_sketch(s, options.kmer_length, options.modulo, {iter_number})))".format(iter_number = i) \
            + ".reduceByKey(lambda a, b: (a + b)).cache()"
        code = code + '\n' + "auxRDD = sketchRDD{iter_number}.filter(lambda v: len(v[1]) >= options.cmax)".format(iter_number = i) \
            + ".flatMap(lambda s: build_aux(s))" \
            + ".reduceByKey(lambda a, b: a + b)"
        code = code + '\n' + "auxList = sc.broadcast(auxRDD.collectAsMap())"
        if(i == 0):
            code = code + '\n' + "outRDD{iter_number} = sketchRDD{iter_number}.filter(lambda v: len(v[1]) > 1 and len(v[1]) < options.cmax)".format(iter_number = i) \
            + ".mapPartitions(break_long_ones)" \
            + ".partitionBy(options.partitions, partitionFunc=rand_part)" \
            + ".mapPartitions(combine_pairs_fast)" \
            + ".aggregateByKey(0, lambda a, b: a + 1, lambda x, y: x + y)" \
            + ".mapPartitions(update_counts).filter(lambda a: a[1] >= options.jmin).persist()"
        else:
            code = code + '\n' + "sketchRDD{iter_number} = sketchRDD{iter_number}.filter(lambda v: len(v[1]) > 1 and len(v[1]) < options.cmax)".format(iter_number = i) \
            + ".mapPartitions(break_long_ones)" \
            + ".partitionBy(options.partitions, partitionFunc=rand_part)" \
            + ".mapPartitions(combine_pairs_fast)" \
            + ".aggregateByKey(0, lambda a, b: a + 1, lambda x, y: x + y)" \
            + ".mapPartitions(update_counts).filter(lambda a: a[1] >= options.jmin)"
            if(i == max(range(options.iterations))):
                code = code + '\n' + "outRDD = outRDD{prev_iter}.union(sketchRDD{iter_number})".format(iter_number = i, prev_iter = i - 1) \
                + ".reduceByKey(lambda a, b: max(a, b)).saveAsTextFile(options.output)"
            else:
                code = code + '\n' + "outRDD{iter_number} = outRDD{prev_iter}.union(sketchRDD{iter_number})".format(iter_number = i, prev_iter = i - 1) \
                + ".reduceByKey(lambda a, b: max(a, b)).persist()" \
                + '\n' + "outRDD{prev_iter}.unpersist()".format(iter_number = i, prev_iter = i - 1)
    with open(options.tmp_file, "w") as f:
        f.write(code)

def compute_similarity_graph(options):
    iteration(options)
    execfile(options.tmp_file)

#Resume GUI section

def initializeOptions():
    '''
    Gather options from command line.
    '''
    print(' ')
    print ("BiSiG Building Tool")
    print ("Copyright (c) 2015-2016 SCoRe Group")
    print(' ')
    usagemessage = "Usage: spark-submit %prog " + usageMessage
    parser = OptionParser(usage=usagemessage)

    #Add parser options
    parser.add_option("-i", "--input", \
        action="store", \
        type="string", \
        dest="ip", \
        default=None, \
        help=inputMessage, \
        metavar="FILE")
    parser.add_option("-o", "--output", \
        action="store", \
        type="string", \
        dest="output", \
        default=None, \
        help=outputMessage, \
        metavar="DIR")
    parser.add_option("-t", "--temp", \
        action="store", \
        type="string", \
        dest="tmp_file", \
        default="tmp.py", \
        help=tmpfileMessage, \
        metavar="FILE")
    parser.add_option("-k", "--kmer", \
        action="store", \
        type="int", \
        dest="kmer_length", \
        default=16, \
        help=kmerlengthMessage, \
        metavar="FILE")
    parser.add_option("-m", "--modulo", \
        action="store", \
        type="int", \
        dest="modulo", \
        default=20, \
        help=moduloMessage)
    parser.add_option("-e", "--iter", \
        action="store", \
        type="int", \
        dest="iterations", \
        default=5, \
        help=iterationsMessage)
    parser.add_option("-c", "--cmax", \
        action="store", \
        type="int", \
        dest="cmax", \
        default=15000, \
        help=cmaxMessage)
    parser.add_option("-j", "--jmin", \
        action="store", \
        type="int", \
        dest="jmin", \
        default=85, \
        help=jminMessage)
    parser.add_option("-p", "--partition", \
        action="store", \
        type="int", \
        dest="partitions", \
        default=1024, \
        help=partitionsMessage)
    parser.add_option("-s", "--smhasher", \
        action="store", \
        type="string", \
        dest="smhasher_loc", \
        default=None, \
        help=smhasherMessage, \
        metavar="DIR")
    (optionvals, sys.argv) = parser.parse_args()
    
    #Checks for...
    #Input
    if (optionvals.ip is None):
        return 2
    try:
        f = open(optionvals.ip, 'r')
        f.close()
    except IOError:
        return 2
    #Output
    if (optionvals.output is None or os.path.isdir(optionvals.output)):
        return 3
    #Kmer Length
    if (optionvals.kmer_length < 2):
        return 4
    #Modulo
    if (optionvals.modulo < 1):
        return 5
    #Iterations
    if (optionvals.iterations < 1):
        return 6
    #Cmax
    if (optionvals.cmax < 1):
        return 7
    #Jmin
    if (optionvals.jmin < 1):
        return 8
    #Number of partitions
    if (optionvals.partitions < 1):
        return 9
    #Smhasher location
    if (optionvals.smhasher_loc is None or not os.path.isdir(optionvals.smhasher_loc + "/lib/python/")):
        return 10
    
    #Assign values to class options.
    options.kmer_length = optionvals.kmer_length
    options.modulo = optionvals.modulo
    options.iterations = optionvals.iterations
    options.cmax = optionvals.cmax
    options.jmin = optionvals.jmin
    options.partitions = optionvals.partitions
    options.tmp_file = optionvals.tmp_file
    options.ip = optionvals.ip
    options.output = optionvals.output
    options.smhasher_loc = optionvals.smhasher_loc
    return 0

#end defs
if (not cmdLineError):
    errcode = initializeOptions()
    if(errcode == 0):
        print('progress can be viewed on most machines at localhost:4040')
        compute_similarity_graph(options)
    else:
        errorMessage(errcode)
