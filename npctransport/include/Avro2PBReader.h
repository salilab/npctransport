/**
 *  \file Avro2PBReader.h
 *  \brief an iterator that reads npctransport protobuf output files one by one
 *         from a list of avro files (used mainly as a wrapper for python usage)
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_AVRO2PB_READER_H
#define IMPNPCTRANSPORT_AVRO2PB_READER_H

#include <avro/ValidSchema.hh>
#include <IMP/npctransport/avro.h>
#include <IMP/npctransport/AvroDataFileData.h>
#include <IMP/base/value_macros.h>
#include <IMP/base/showable_macros.h>
#include <avro/DataFile.hh>
#include <fstream>
#include <iomanip>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT Avro2PBReader {
 private:
  typedef avro::DataFileReader<IMP_npctransport::wrapper> t_avro_reader;

 public:
  /** Initiates a reader that goes over all output
      entries in all files specified in avro_filenames
  */
  Avro2PBReader(std::vector<std::string> avro_filenames);

  /** Initiates a reader that goes over all output
      entries in all files specified in avro_filenames
  */
  Avro2PBReader(std::string avro_filename);


  /** closes any open files */
  ~Avro2PBReader();

  /**
     Read the next output entry into output and returns it
     as string. If no input is left, returns "" and invalidates
     this object.
  */
  std::string read_next();

  //! returns true if there are still files to go over
  //! (though possibly no entries left in neither of them)
  bool get_is_valid();

 private:
  //! close any open file if one exists and move cursor to next file index
  void advance_current_reader();

  // called from ctr, this pain is needed since constructor delegation
  // is only supported from g++ 4.7, so we use init() for backward compatibility
  void init(std::vector<std::string> avro_filenames);

 private:
  std::vector<std::string> avro_filenames_; // list of files to go over
  t_avro_reader* avro_reader_;
  unsigned int cur_file_; // file index we're reading now

 public:
  IMP_SHOWABLE_INLINE(Avro2PBReader,
                      out << "Avro2PBReader with "
                           << avro_filenames_.size() << " input avro files"  );

};

IMP_VALUES(Avro2PBReader, Avro2PBReaders);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_AVRO2PB_READER_H */
