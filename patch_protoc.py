"""Patch the protoc-generated header to resolve export symbols"""

import sys
import os
import shutil

header_fn = sys.argv[1]
link_header_fn = sys.argv[2]

if sys.version_info[0] == 2:
    linesep = os.linesep
else:
    linesep = bytes(os.linesep, 'ascii')

with open(header_fn, 'rb') as fh:
    contents = fh.read()

with open(header_fn, 'wb') as fh:
    fh.write(b'#include <IMP/npctransport/npctransport_config.h>' + linesep)
    # Don't warn about missing override keyword; this is not our code
    fh.write(b'#if __GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 1)'
             + linesep)
    fh.write(b'IMP_GCC_PRAGMA(diagnostic ignored "-Wsuggest-override")'
             + linesep)
    fh.write(b'#endif' + linesep)
    fh.write(contents)

if os.path.exists(link_header_fn):
    os.unlink(link_header_fn)

# symlinks usually don't work on Windows, so copy instead
if sys.platform == 'win32':
    shutil.copy(header_fn, link_header_fn)
else:
    os.symlink(header_fn, link_header_fn)
