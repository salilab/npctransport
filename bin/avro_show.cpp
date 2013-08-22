/**
 * \file avro_cat.cpp
 * \brief Concatenate a bunch of protobufs to an avro archive
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#include <IMP/npctransport/avro.h>
#include <IMP/npctransport/AvroDataFileData.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <DataFile.hh>
#include <google/protobuf/text_format.h>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " avrofiles" << std::endl;
    return 1;
  }
  std::string out = argv[argc - 1];

  for (int i = 1; i < argc; ++i) {
    avro::DataFileReader<IMP_npctransport::wrapper> rd(
        argv[i], IMP::npctransport::get_avro_data_file_schema());
    rd.readerSchema().toJson(std::cout);
    rd.dataSchema().toJson(std::cout);
    IMP_npctransport::wrapper data;
    while (rd.read(data)) {
      std::string str(data.value.begin(), data.value.end());
      npctransport_proto::Output output;
      output.ParseFromString(str);
      std::string output_string;
      google::protobuf::TextFormat::PrintToString(output, &output_string);
      std::cout << "work unit: " << data.key << std::endl;
      std::cout << output_string << std::endl;
    }
  }
  return 0;
}
