/**
 *  \file Avro2PBReader.cc
 *  \brief an iterator that reads npctransport protobuf output files one by one
 *         from a list of avro files (used mainly as a wrapper for python usage)
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/Avro2PBReader.h>
#include <ValidSchema.hh>
#include <IMP/npctransport/avro.h>
#include <IMP/npctransport/AvroDataFileData.h>
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
IMP_GCC_PRAGMA(diagnostic push)
#endif
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "npctransport.pb.h"
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
IMP_GCC_PRAGMA(diagnostic pop)
#endif

#include <DataFile.hh>
//#include <google/protobuf/text_format.h>
#include <IMP/exception.h>
#include <fstream>
#include <iomanip>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

Avro2PBReader::Avro2PBReader(const Strings& avro_filenames) {
  init(avro_filenames);
}

Avro2PBReader::Avro2PBReader(std::string avro_filename)
    // TODO: from gcc 4.7 it is possible to delegate
    // constructors but for now we keep backward comapibility pain
    {
  init(Strings(1, avro_filename));
}

void Avro2PBReader::init(const Strings& avro_filenames) {
  avro_filenames_ = avro_filenames;
  avro_reader_ = nullptr;
  cur_file_ = 0;
}

/** closes any open files */
Avro2PBReader::~Avro2PBReader() { advance_current_reader(); }

std::string Avro2PBReader::read_next() {
  if (!get_is_valid()) {
    return "";
  }
  if (!avro_reader_) {
    avro_reader_ =
        new t_avro_reader(avro_filenames_[cur_file_].c_str(),
                          IMP::npctransport::get_avro_data_file_schema());
  }
  IMP_npctransport::wrapper data;

  if (!(avro_reader_->read(data))) {  // no more data = go to next
    advance_current_reader();
    return read_next();
  }
  return std::string(data.value.begin(), data.value.end());
}

bool Avro2PBReader::get_is_valid() {
  bool is_valid = (cur_file_ < avro_filenames_.size());
  return is_valid;
}

/*************** Private ****************/

void Avro2PBReader::advance_current_reader() {
  if (avro_reader_) delete avro_reader_;
  avro_reader_ = nullptr;
  cur_file_++;
}

IMPNPCTRANSPORT_END_NAMESPACE
