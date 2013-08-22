/**
 * \file avro_cat.cpp
 * \brief Concatenate a bunch of protobufs to an avro archive
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#include <IMP/npctransport/avro.h>
#include <IMP/npctransport/AvroDataFileData.h>
#include <DataFile.hh>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " protobufs... avrofile" << std::endl;
    return 1;
  }
  std::string out = argv[argc - 1];

  avro::DataFileWriter<IMP_npctransport::wrapper> wr(
      out.c_str(), IMP::npctransport::get_avro_data_file_schema());

  for (int i = 0; i < argc - 2; ++i) {
    IMP_npctransport::wrapper data;
    data.key = "none";
    std::ifstream file(argv[i],
                       std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
      throw std::runtime_error("couldn't open file");
    }
    data.value.resize(file.tellg());

    file.seekg(0, std::ios::beg);
    if (!file.read(reinterpret_cast<char *>(&data.value[0]), data.value.size()))
      wr.write(data);
  }
  return 0;
}
