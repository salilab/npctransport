message(STATUS "Setting up proto " "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.cpp")

# there is a #include 'npctransport.ph.h' in the cpp file
add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.cpp"
                          COMMAND protoc "--cpp_out=dllexport_decl=IMPNPCTRANSPORTEXPORT:."
                          "-I${CMAKE_SOURCE_DIR}/modules/npctransport/data/"
                          "${CMAKE_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMAND mv npctransport.pb.cc npctransport.pb.cpp
                          # add config header to resolve export symbols
                          COMMAND mv "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h" "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h.out"
                          COMMAND echo "\"\#include\"" "\"<IMP/npctransport/npctransport_config.h>\"" > "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          COMMAND cat "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h.out" >> "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          COMMAND rm -f "${CMAKE_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          COMMAND ln -s "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h" "${CMAKE_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          DEPENDS "${CMAKE_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMENT "Creating protoc stuff for npctransport"
                          WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/npctransport")

add_custom_target(IMP.npctransport-proto ALL DEPENDS
                          "${CMAKE_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.cpp")
set_property(TARGET IMP.npctransport-proto PROPERTY FOLDER "IMP.npctransport")

list(APPEND IMP_NPCTRANSPORT_LIBRARY_EXTRA_SOURCES "${CMAKE_BINARY_DIR}/src/npctransport/npctransport.pb.cpp")

set(IMP_NPCTRANSPORT_LIBRARY_EXTRA_DEPENDENCIES "IMP.npctransport-proto" CACHE INTERNAL "" FORCE)


# there is a #include 'npctransport.ph.h' in the cpp file
add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/lib/IMP/npctransport/npctransport_pb2.py"
                          COMMAND protoc "--python_out=."
                          "-I${CMAKE_SOURCE_DIR}/modules/npctransport/data/"
                          "${CMAKE_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMAND mv npctransport_pb2.py "${CMAKE_BINARY_DIR}/lib/IMP/npctransport/"
                          # add config header to resolve export symbols
                          DEPENDS "${CMAKE_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMENT "Creating python protoc stuff for npctransport"
                          WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/npctransport")

add_custom_target(IMP.npctransport-python_proto ALL DEPENDS "${CMAKE_BINARY_DIR}/lib/IMP/npctransport/npctransport_pb2.py" )
set_property(TARGET IMP.npctransport-python_proto PROPERTY FOLDER "IMP.npctransport")

include_directories(${CMAKE_BINARY_DIR}/include/IMP/npctransport/internal)

set(IMP_NPCTRANSPORT_CONFIG "IMP_NPCTRANSPORT_AVRO_NAMESPACE=rmf_avro:IMP_NPCTRANSPORT_PROTOBUF_NAMESPACE=google\\:\\:protobuf")

set(IMP_NPCTRANSPORT_PYTHON_EXTRA_DEPENDENCIES IMP.npctransport-python_proto CACHE INTERNAL "" FORCE)
