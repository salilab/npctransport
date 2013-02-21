message(STATUS "Setting up avro "${PROJECT_BINARY_DIR}/include/IMP/npctransport/AvroDataFileData.h)

add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/include/IMP/npctransport/AvroDataFileData.h"
     COMMAND avrogencpp "--input" "${PROJECT_SOURCE_DIR}/modules/npctransport/data/AvroDataFileData.json"
     "--output" "${PROJECT_BINARY_DIR}/include/IMP/npctransport/AvroDataFileData.h" "--namespace" "IMP_npctransport"
     DEPENDS "${PROJECT_SOURCE_DIR}/modules/npctransport/data/AvroDataFileData.json" COMMENT "Creating json header for npctransport")

add_custom_target(npctransport_avro ALL DEPENDS "${PROJECT_BINARY_DIR}/include/IMP/npctransport/AvroDataFileData.h" )

message(STATUS "Setting up proto " "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.cpp")

# there is a #include 'npctransport.ph.h' in the cpp file
add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/include/IMP/npctransport/internal/npctransport.pb.h"
                          "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.h"
                          "${PROJECT_BINARY_DIR}/src/npctransport/npctransport.pb.cpp"
                          COMMAND mkdir -p "${PROJECT_BINARY_DIR}/src/npctransport/"
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