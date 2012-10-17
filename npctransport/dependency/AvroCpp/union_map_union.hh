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


#ifndef UNION_MAP_UNION_HH_4104945341__H_
#define UNION_MAP_UNION_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace umu {
struct _union_map_union_Union__0__ {
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
    _union_map_union_Union__0__();
};

struct r3 {
    std::string name;
    std::vector<uint8_t> data;
    _union_map_union_Union__0__ rev;
};

struct _union_map_union_Union__1__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    std::map<std::string, r3 > get_map() const;
    void set_map(const std::map<std::string, r3 >& v);
    bool is_null() const {
        return (idx_ == 1);
    }
    void set_null() {
        idx_ = 1;
        value_ = boost::any();
    }
    _union_map_union_Union__1__();
};

struct r1 {
    std::string id;
    _union_map_union_Union__1__ val;
};

inline
std::string _union_map_union_Union__0__::get_string() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::string >(value_);
}

inline
void _union_map_union_Union__0__::set_string(const std::string& v) {
    idx_ = 0;
    value_ = v;
}

inline
std::map<std::string, r3 > _union_map_union_Union__1__::get_map() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::map<std::string, r3 > >(value_);
}

inline
void _union_map_union_Union__1__::set_map(const std::map<std::string, r3 >& v) {
    idx_ = 0;
    value_ = v;
}

inline _union_map_union_Union__0__::_union_map_union_Union__0__() : idx_(0), value_(std::string()) { }
inline _union_map_union_Union__1__::_union_map_union_Union__1__() : idx_(0), value_(std::map<std::string, r3 >()) { }
}
namespace avro {
template<> struct codec_traits<umu::_union_map_union_Union__0__> {
    static void encode(Encoder& e, umu::_union_map_union_Union__0__ v) {
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
    static void decode(Decoder& d, umu::_union_map_union_Union__0__& v) {
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

template<> struct codec_traits<umu::r3> {
    static void encode(Encoder& e, const umu::r3& v) {
        avro::encode(e, v.name);
        avro::encode(e, v.data);
        avro::encode(e, v.rev);
    }
    static void decode(Decoder& d, umu::r3& v) {
        avro::decode(d, v.name);
        avro::decode(d, v.data);
        avro::decode(d, v.rev);
    }
};

template<> struct codec_traits<umu::_union_map_union_Union__1__> {
    static void encode(Encoder& e, umu::_union_map_union_Union__1__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            avro::encode(e, v.get_map());
            break;
        case 1:
            e.encodeNull();
            break;
        }
    }
    static void decode(Decoder& d, umu::_union_map_union_Union__1__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            {
                std::map<std::string, umu::r3 > vv;
                avro::decode(d, vv);
                v.set_map(vv);
            }
            break;
        case 1:
            d.decodeNull();
            v.set_null();
            break;
        }
    }
};

template<> struct codec_traits<umu::r1> {
    static void encode(Encoder& e, const umu::r1& v) {
        avro::encode(e, v.id);
        avro::encode(e, v.val);
    }
    static void decode(Decoder& d, umu::r1& v) {
        avro::decode(d, v.id);
        avro::decode(d, v.val);
    }
};

}
#endif
