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


#ifndef UNION_CONFLICT_HH_4104945341__H_
#define UNION_CONFLICT_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace uc {
struct _union_conflict_Union__0__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    std::string get_string() const;
    void set_string(const std::string& v);
    bool is_null() const {
        return (idx_ == 1);
    }
    void set_null() {
        idx_ = 1;
        value_ = boost::any();
    }
    _union_conflict_Union__0__();
};

struct uc {
    std::string rev_t;
    std::vector<uint8_t> data;
    _union_conflict_Union__0__ rev;
};

inline
std::string _union_conflict_Union__0__::get_string() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::string >(value_);
}

inline
void _union_conflict_Union__0__::set_string(const std::string& v) {
    idx_ = 0;
    value_ = v;
}

inline _union_conflict_Union__0__::_union_conflict_Union__0__() : idx_(0), value_(std::string()) { }
}
namespace avro {
template<> struct codec_traits<uc::_union_conflict_Union__0__> {
    static void encode(Encoder& e, uc::_union_conflict_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            avro::encode(e, v.get_string());
            break;
        case 1:
            e.encodeNull();
            break;
        }
    }
    static void decode(Decoder& d, uc::_union_conflict_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            {
                std::string vv;
                avro::decode(d, vv);
                v.set_string(vv);
            }
            break;
        case 1:
            d.decodeNull();
            v.set_null();
            break;
        }
    }
};

template<> struct codec_traits<uc::uc> {
    static void encode(Encoder& e, const uc::uc& v) {
        avro::encode(e, v.rev_t);
        avro::encode(e, v.data);
        avro::encode(e, v.rev);
    }
    static void decode(Decoder& d, uc::uc& v) {
        avro::decode(d, v.rev_t);
        avro::decode(d, v.data);
        avro::decode(d, v.rev);
    }
};

}
#endif
