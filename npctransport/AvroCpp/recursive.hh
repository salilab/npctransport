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


#ifndef RECURSIVE_HH_4104945341__H_
#define RECURSIVE_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace rec {
struct LongList;
struct _recursive_Union__0__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    LongList get_LongList() const;
    void set_LongList(const LongList& v);
    bool is_null() const {
        return (idx_ == 1);
    }
    void set_null() {
        idx_ = 1;
        value_ = boost::any();
    }
    _recursive_Union__0__();
};

struct LongList {
    int64_t value;
    _recursive_Union__0__ next;
};

inline
LongList _recursive_Union__0__::get_LongList() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<LongList >(value_);
}

inline
void _recursive_Union__0__::set_LongList(const LongList& v) {
    idx_ = 0;
    value_ = v;
}

inline _recursive_Union__0__::_recursive_Union__0__() : idx_(0), value_(LongList()) { }
}
namespace avro {
template<> struct codec_traits<rec::_recursive_Union__0__> {
    static void encode(Encoder& e, rec::_recursive_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            avro::encode(e, v.get_LongList());
            break;
        case 1:
            e.encodeNull();
            break;
        }
    }
    static void decode(Decoder& d, rec::_recursive_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            {
                rec::LongList vv;
                avro::decode(d, vv);
                v.set_LongList(vv);
            }
            break;
        case 1:
            d.decodeNull();
            v.set_null();
            break;
        }
    }
};

template<> struct codec_traits<rec::LongList> {
    static void encode(Encoder& e, const rec::LongList& v) {
        avro::encode(e, v.value);
        avro::encode(e, v.next);
    }
    static void decode(Decoder& d, rec::LongList& v) {
        avro::decode(d, v.value);
        avro::decode(d, v.next);
    }
};

}
#endif
