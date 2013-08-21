#!/usr/bin/python

import sys
import column_num

def usage_and_exit():
    print "Usage:"
    print sys.argv[0], " <filename> <caption1> [<caption2> ...]"
    print " ... |", sys.argv[0], "<caption1> [<caption2> ...]"
    print
    print "Prints the columns with specified captions in the header line of"
    print  " <filename>, or stdin if unix piping"
    exit(-1)

# getter for multiple list elements
getListElements = lambda searchList, ind: [searchList[i] for i in ind]

i= 1
# open file or stdin pipe as file
file= None
if(len(sys.argv)<=1): usage_and_exit()
if(not sys.stdin.isatty()): # stdin not from terminal = piping
    file= sys.stdin
else: # else use file input
    if(len(sys.argv)<=2): usage_and_exit()
    filename= sys.argv[i]
    file=open(filename, 'r')
    i=i+1
# get indexes of captions
captions= sys.argv[i:]
header_list = file.readline().split()
indexes = column_num.get_column_indexes(header_list, captions)
# print captions and data
print indexes
for caption in captions: print caption,
print
for line in file.readlines():
    line_split = line.split()
    for num in getListElements(line_split, indexes):
        print num,
    print
file.close()
