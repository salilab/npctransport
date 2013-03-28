message(STATUS "Setting up proto " "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.cpp")

# there is a #include 'npctransport.ph.h' in the cpp file
add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.cpp"
                          COMMAND protoc "--cpp_out=dllexport_decl=IMPNPCTRANSPORTEXPORT:."
                          "-I${PROJECT_SOURCE_DIR}/modules/npctransport/data/"
                          "${PROJECT_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMAND mv npctransport.pb.cc npctransport.pb.cpp
                          # add config header to resolve export symbols
                          COMMAND mv "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h" "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h.out"
                          COMMAND echo "\"\#include\"" "\"<IMP/npctransport/npctransport_config.h>\"" > "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          COMMAND cat "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h.out" >> "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          COMMAND ln -s "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h" "${PROJECT_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          DEPENDS "${PROJECT_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMENT "Creating protoc stuff for npctransport"
                          WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src/npctransport")

add_custom_target(npctransport_proto ALL DEPENDS "${PROJECT_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h" )

set(IMP_NPCTRANSPORT_LIBRARY_EXTRA_SOURCES "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.cpp" "${PROJECT_BINARY_DIR}/include/IMP/npctransport/AvroDataFileData.h" CACHE INTERNAL "" FORCE)

set(IMP_NPCTRANSPORT_LIBRARY_EXTRA_DEPENDENCIES npctransport_proto npctransport_avro CACHE INTERNAL "" FORCE)


# there is a #include 'npctransport.ph.h' in the cpp file
add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/lib/IMP/npctransport/npctransport_pb2.py"
                          COMMAND protoc "--python_out=."
                          "-I${PROJECT_SOURCE_DIR}/modules/npctransport/data/"
                          "${PROJECT_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMAND mv npctransport_pb2.py "${PROJECT_BINARY_DIR}/lib/IMP/npctransport/"
                          # add config header to resolve export symbols
                          DEPENDS "${PROJECT_SOURCE_DIR}/modules/npctransport/data/npctransport.proto"
                          COMMENT "Creating python protoc stuff for npctransport"
                          WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src/npctransport")

add_custom_target(npctransport_python_proto ALL DEPENDS "${PROJECT_BINARY_DIR}/lib/IMP/npctransport/npctransport_pb2.py" )

set(IMP_NPCTRANSPORT_PYTHON_EXTRA_DEPENDENCIES npctransport_proto CACHE INTERNAL "" FORCE)
