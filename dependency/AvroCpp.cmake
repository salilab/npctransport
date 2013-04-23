
message(STATUS "Using internal AvroCPP")
set(AVROCPP_BINARY_DIR ${PROJECT_BINARY_DIR}/src/dependency/AvroCpp CACHE INTERNAL "" FORCE)
add_subdirectory("${PROJECT_SOURCE_DIR}/modules/npctransport/dependency/AvroCpp" "${AVROCPP_BINARY_DIR}")
set(AVROCPP_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/modules/npctransport/dependency/AvroCpp/" CACHE INTERNAL "" FORCE)
set(AVROCPP_LIBRARIES avrocpp CACHE INTERNAL "" FORCE)
file(WRITE "${PROJECT_BINARY_DIR}/data/build_info/AvroCpp" "ok=True
includepath=\"${AVROCPP_INCLUDE_PATH}\"
")
