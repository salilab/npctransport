
if (${RMF_AVRO} MATCHES "internal")
  message(STATUS "Using internal AvroCPP")
  set(AVROGENCPP ${RMF_BINARY_DIR}/AvroCpp/avrogencpp CACHE INTERNAL "" FORCE)
  set(AVROGENCPP_DEPENDENCY ${RMF_BINARY_DIR}/AvroCpp/avrogencpp CACHE INTERNAL "" FORCE)
  set(AVROCPP_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/modules/rmf/dependency/AvroCpp CACHE INTERNAL "" FORCE)
  set(ACROCPP_LIBRARIES avrocpp CACHE INTERNAL "" FORCE)
  file(WRITE "${PROJECT_BINARY_DIR}/data/build_info/AvroCpp" "ok=True
includepath=\"${AvroCpp_INCLUDE_PATH}\"
")
else()
  set(AVROGENCPP avrogencpp CACHE INTERNAL "" FORCE)
endif()