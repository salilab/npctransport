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


#ifndef CIRCULARDEP_HH_4104945341__H_
#define CIRCULARDEP_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace cd {
struct Information;
struct _circulardep_Union__0__ {
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
    Information get_Information() const;
    void set_Information(const Information& v);
    _circulardep_Union__0__();
};

struct Item {
    std::string id;
    _circulardep_Union__0__ entities;
};

struct _circulardep_Union__1__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    int32_t get_int() const;
    void set_int(const int32_t& v);
    double get_double() const;
    void set_double(const double& v);
    _circulardep_Union__1__();
};

struct Information {
    std::string id;
    Item externalItem;
    _circulardep_Union__1__ innerUnion;
};

inline
Information _circulardep_Union__0__::get_Information() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<Information >(value_);
}

inline
void _circulardep_Union__0__::set_Information(const Information& v) {
    idx_ = 1;
    value_ = v;
}

inline
int32_t _circulardep_Union__1__::get_int() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<int32_t >(value_);
}

inline
void _circulardep_Union__1__::set_int(const int32_t& v) {
    idx_ = 0;
    value_ = v;
}

inline
double _circulardep_Union__1__::get_double() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<double >(value_);
}

inline
void _circulardep_Union__1__::set_double(const double& v) {
    idx_ = 1;
    value_ = v;
}

inline _circulardep_Union__0__::_circulardep_Union__0__() : idx_(0) { }
inline _circulardep_Union__1__::_circulardep_Union__1__() : idx_(0), value_(int32_t()) { }
}
namespace avro {
template<> struct codec_traits<cd::_circulardep_Union__1__> {
    static void encode(Encoder& e, cd::_circulardep_Union__1__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            avro::encode(e, v.get_int());
            break;
        case 1:
            avro::encode(e, v.get_double());
            break;
        }
    }
    static void decode(Decoder& d, cd::_circulardep_Union__1__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            {
                int32_t vv;
                avro::decode(d, vv);
                v.set_int(vv);
            }
            break;
        case 1:
            {
                double vv;
                avro::decode(d, vv);
                v.set_double(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<cd::Information> {
    static void encode(Encoder& e, const cd::Information& v) {
        avro::encode(e, v.id);
        avro::encode(e, v.externalItem);
        avro::encode(e, v.innerUnion);
    }
    static void decode(Decoder& d, cd::Information& v) {
        avro::decode(d, v.id);
        avro::decode(d, v.externalItem);
        avro::decode(d, v.innerUnion);
    }
};

template<> struct codec_traits<cd::_circulardep_Union__0__> {
    static void encode(Encoder& e, cd::_circulardep_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_Information());
            break;
        }
    }
    static void decode(Decoder& d, cd::_circulardep_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                cd::Information vv;
                avro::decode(d, vv);
                v.set_Information(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<cd::Item> {
    static void encode(Encoder& e, const cd::Item& v) {
        avro::encode(e, v.id);
        avro::encode(e, v.entities);
    }
    static void decode(Decoder& d, cd::Item& v) {
        avro::decode(d, v.id);
        avro::decode(d, v.entities);
    }
};

}
#endif
