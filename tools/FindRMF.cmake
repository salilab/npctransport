#[=======================================================================[.rst:
FindRMF
-------

Try to find RMF. This should be called after finding IMP and *before*
including $IMP_USE_FILE, as RMF is often bundled with IMP (and FindIMP
will overwrite RMF variables).

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``RMF_FOUND``
  system has RMF
``RMF_INCLUDE_PATH``
  the RMF include directory
``RMF_SWIG_DIR``
  the directory containing SWIG (.i) files for RMF
``RMF_LIBRARIES``
  Link this to use RMF
``RMF_VERSION_STRING``
  the version of RMF found


#]=======================================================================]

# If paths were already provided by IMP but are wrong, clear them so we
# can try again
if (NOT EXISTS ${RMF_SWIG_DIR}/RMF.i)
  unset(RMF_SWIG_DIR)
  unset(RMF_LIBRARIES)
  unset(RMF_INCLUDE_PATH)
endif()

find_path(RMF_INCLUDE_PATH RMF.h PATH_SUFFIXES include)
find_path(RMF_SWIG_DIR RMF.i PATH_SUFFIXES share/RMF/swig)
if (NOT RMF_LIBRARIES)
  find_library(RMF_LIBRARIES NAMES RMF PATH_SUFFIXES lib)
endif()

if (RMF_INCLUDE_PATH AND EXISTS "${RMF_INCLUDE_PATH}/RMF/config.h")
  file(STRINGS "${RMF_INCLUDE_PATH}/RMF/config.h" RMF_MAJOR_H REGEX "#define RMF_VERSION_MAJOR +([0-9]+)")
  file(STRINGS "${RMF_INCLUDE_PATH}/RMF/config.h" RMF_MINOR_H REGEX "#define RMF_VERSION_MINOR +([0-9]+)")
  string(REGEX REPLACE " *#define RMF_VERSION_MAJOR +([0-9]+) *" "\\1" RMF_VERSION_MAJOR "${RMF_MAJOR_H}")
  string(REGEX REPLACE " *#define RMF_VERSION_MINOR +([0-9]+) *" "\\1" RMF_VERSION_MINOR "${RMF_MINOR_H}")
  set(RMF_VERSION_STRING "${RMF_VERSION_MAJOR}.${RMF_VERSION_MINOR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RMF
        REQUIRED_VARS RMF_LIBRARIES RMF_INCLUDE_PATH RMF_SWIG_DIR
        VERSION_VAR RMF_VERSION_STRING)
