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


#ifndef UNION_ARRAY_UNION_HH_4104945341__H_
#define UNION_ARRAY_UNION_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace uau {
struct _union_array_union_Union__0__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    bool is_null() const {
        return (idx_ == 0);
    }
    void set_null() {
        idx_ = 0;
        value_ = boost::any();
    }
    int32_t get_int() const;
    void set_int(const int32_t& v);
    _union_array_union_Union__0__();
};

struct _union_array_union_Union__1__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    bool is_null() const {
        return (idx_ == 0);
    }
    void set_null() {
        idx_ = 0;
        value_ = boost::any();
    }
    std::vector<_union_array_union_Union__0__ > get_array() const;
    void set_array(const std::vector<_union_array_union_Union__0__ >& v);
    _union_array_union_Union__1__();
};

struct r1 {
    _union_array_union_Union__1__ f1;
};

inline
int32_t _union_array_union_Union__0__::get_int() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<int32_t >(value_);
}

inline
void _union_array_union_Union__0__::set_int(const int32_t& v) {
    idx_ = 1;
    value_ = v;
}

inline
std::vector<_union_array_union_Union__0__ > _union_array_union_Union__1__::get_array() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::vector<_union_array_union_Union__0__ > >(value_);
}

inline
void _union_array_union_Union__1__::set_array(const std::vector<_union_array_union_Union__0__ >& v) {
    idx_ = 1;
    value_ = v;
}

inline _union_array_union_Union__0__::_union_array_union_Union__0__() : idx_(0) { }
inline _union_array_union_Union__1__::_union_array_union_Union__1__() : idx_(0) { }
}
namespace avro {
template<> struct codec_traits<uau::_union_array_union_Union__0__> {
    static void encode(Encoder& e, uau::_union_array_union_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_int());
            break;
        }
    }
    static void decode(Decoder& d, uau::_union_array_union_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                int32_t vv;
                avro::decode(d, vv);
                v.set_int(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<uau::_union_array_union_Union__1__> {
    static void encode(Encoder& e, uau::_union_array_union_Union__1__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_array());
            break;
        }
    }
    static void decode(Decoder& d, uau::_union_array_union_Union__1__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                std::vector<uau::_union_array_union_Union__0__ > vv;
                avro::decode(d, vv);
                v.set_array(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<uau::r1> {
    static void encode(Encoder& e, const uau::r1& v) {
        avro::encode(e, v.f1);
    }
    static void decode(Decoder& d, uau::r1& v) {
        avro::decode(d, v.f1);
    }
};

}
#endif
