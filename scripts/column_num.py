#!/usr/bin/python

import sys

def usage_and_exit():
    print "Usage:"
    print sys.argv[0], " <filename> <caption1> [<caption2> ...]"
    print " ... |", sys.argv[0], "<caption1> [<caption2> ...]"
    print
    print "Prints the indexes of specified captions in the header line of"
    print  " <filename>, or stdin if unix piping"
    exit(-1)

def get_column_indexes(header_list, captions):
    ret=[]
    for caption in captions:
        try:
            k=header_list.index(caption)
            ret.append(k)
        except ValueError:
            ret.append('NA')
    return ret

if __name__ == "__main__":
    i= 1
    file= None
    if(len(sys.argv)<=1): usage_and_exit()
    if(not sys.stdin.isatty()): # stdin not from terminal = piping
        file= sys.stdin
    else: # else use file input
        if(len(sys.argv)<=2): usage_and_exit()
        filename= sys.argv[i]
        file=open(filename, 'r')
        i=i+1
    # the rest of the arguments are the column captions
    captions= sys.argv[i:]
    header_list = file.readline().split()
    indexes = get_column_indexes(header_list, captions)
    for index in indexes: print index,
    print
    file.close()
