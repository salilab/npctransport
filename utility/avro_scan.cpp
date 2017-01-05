/**
 * \file avro_cat.cpp
 * \brief Concatenate a bunch of protobufs to an avro archive
 * Copyright 2007-2017 IMP Inventors. All rights reserved.
 */
#include <IMP/npctransport/avro.h>
#include <IMP/npctransport/AvroDataFileData.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/check_macros.h>
#include <DataFile.hh>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " avrofile" << std::endl;
    return 1;
  }
  std::string out = argv[argc - 1];

  IMP_NPCTRANSPORT_AVRO_NAMESPACE::DataFileReader<IMP_npctransport::wrapper> rd(
      out.c_str(), IMP::npctransport::get_avro_data_file_schema());

  IMP_npctransport::wrapper data;
  while (rd.read(data)) {
    std::cout << data.key << " ";
    npctransport_proto::Output pb;
    std::string str(data.value.begin(), data.value.end());
    IMP_INTERNAL_CHECK(
        str.size() == data.value.size(),
        "Sizes don't match " << str.size() << " vs " << data.value.size());
    pb.ParseFromString(str);
    std::cout << " work unit " << pb.assignment().work_unit() << std::endl;
  }
  return 0;
}
