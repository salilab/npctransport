set(AVROCPP_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/modules/rmf/dependency/RMF/src/avrocpp/api CACHE INTERNAL "" FORCE)
set(AVROCPP_LIBRARIES "" CACHE INTERNAL "" FORCE)
list(APPEND IMP_npctransport_LIBRARY_EXTRA_SOURCES avrocpp)
file(WRITE "${CMAKE_BINARY_DIR}/data/build_info/AvroCpp" "ok=True\nincludepath=\"${ACROCPP_INCLUDE_PATH}\"\nlibpath=\"${AVROCPP_LIBRARIES}\"\n")
