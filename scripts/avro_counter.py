# #!/netapp/sali/barak/Runs/NPC_BD/Cyl2/tools/imppy.sh /usr/bin/python
from IMP.npctransport import *
import numpy
import math
import sys
import re
import numbers
import os
import glob
import traceback
KEY_CAPTIONS = None

def query_yes_no(question, default=None):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")



def get_file_paths_generator(files_and_folders):
    """
    return a generator that iterates over:
    1) all the files that are not folders in files_and_folders
    2) all the inner *.pb files for every folder in files_and_folders
    """
    for path in files_and_folders:
        if(not os.path.exists(path)):
            print "WARNING: input path", path, "does not exist"
            continue
        if(os.path.isdir(path)):
            print "DIR PATH", path
            pb_paths = glob.glob(path + "/*.pb")
            for pb_path in pb_paths:
                print "PB_PATH", pb_path
                yield pb_path
        else:
            print "FILE_PATH", path
            yield path

def get_pb_messages_generator(files_and_folders, parse = False):
    """
    return a generator that iterates over all pb output messages in:
    1) all the non-avro and avro files in files_and_folders
    2) all the inner *.pb for every folder in files_and_folders

    The generator returns a sequence of tuples (message, tag)
    with message the protobuf message if parse=True or the raw
    message string if parse==False, and the tag is the file name
    where the message is stored, with an internal id for avro files
    """
    o = Output()
    for file_path in get_file_paths_generator(files_and_folders):
        prefix,ext=os.path.splitext(file_path)
        if ext==".pb":
            try:
                FILE=open(file_path,"rb")
                if(parse):
                    o.ParseFromString(FILE.read())
                    yield o, file_path
                else:
                    yield FILE.read(), file_path
            except KeyboardInterrupt: exit(-1)
            except:
                print >> sys.stderr, 'Unexpected error: file %s' % file_name, sys.exc_info()
                continue
        else: # assume avro
            avro_reader=Avro2PBReader([file_path])
            i=0
            nerr=0
            MAX_ERRORS_PER_AVRO = 50
            while(avro_reader.get_is_valid()):
                #  print a.is_valid()
                try:
                    i=i+1
                    s = avro_reader.read_next()
                    if(s == ""): break; # invalid output = the end
                    tag = file_path + "." + str(i)
                    if(parse):
                        o.ParseFromString(s)
                        yield o, tag
                    else:
                        yield s, tag
                except KeyboardInterrupt: exit(-1)
                except:
                    nerr = nerr + 1
                    print >> sys.stderr, 'Unexpected error: file %s error %d' \
                        % (avro_reader.get_cur_file_name(), nerr), sys.exc_info()
                    if(nerr <= MAX_ERRORS_PER_AVRO):
                        print >> sys.stderr, "CONTINUING"
                        continue
                    else:
                        print >> sys.stderr, "BREAKING from file %s due to excessive errors" \
                            % avro_reader.get_cur_file_name()
                        break






#################################
'''
Read pb or avro files in list of files and folders sys.arv[2:],
make stats on them, and write them to sys.argb[1]
'''
if __name__ != "__main__": exit()
files_and_folders = sys.argv[1:]
print "Parsing %d files and folders" % len(files_and_folders)
results = {}
n_total = 0
k_flush=1000
is_parse = False
for pb_message, tag in get_pb_messages_generator(files_and_folders, is_parse):
    try:
        n_total = n_total + 1
        if(n_total % k_flush == 0):
            print "========== n = %d =========" % n_total
    except:
        type,msg,tb = sys.exc_info()
        print >> sys.stderr, 'Unexpected error in message tag %s - %s %s' \
            % (tag, type, msg)
        traceback.print_tb(tb)
# FINAL OUTPUT
print "========== n = %d =========" % n_total
