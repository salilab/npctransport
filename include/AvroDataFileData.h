/**
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef INCLUDE_ACRODATAFILEDATA_H_2861103865__H_
#define INCLUDE_ACRODATAFILEDATA_H_2861103865__H_

#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace IMP_npctransport {
struct wrapper {
  std::string key;
  std::vector<uint8_t> value;
};
}
namespace IMP_NPCTRANSPORT_AVRO_NAMESPACE {
template <>
struct codec_traits<IMP_npctransport::wrapper> {
  template <class Encoder>
  static void encode(Encoder& e, const IMP_npctransport::wrapper& v) {
    IMP_NPCTRANSPORT_AVRO_NAMESPACE::encode(e, v.key);
    IMP_NPCTRANSPORT_AVRO_NAMESPACE::encode(e, v.value);
  }
  template <class Decoder>
  static void decode(Decoder& d, IMP_npctransport::wrapper& v) {
    IMP_NPCTRANSPORT_AVRO_NAMESPACE::decode(d, v.key);
    IMP_NPCTRANSPORT_AVRO_NAMESPACE::decode(d, v.value);
  }
};
}
#endif
