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


#ifndef BIGRECORD2_HH_4104945341__H_
#define BIGRECORD2_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace testgen2 {
struct Nested {
    int32_t inval3;
    std::string inval2;
    double inval1;
};

enum ExampleEnum {
    three,
    two,
    one,
    zero,
};

struct _bigrecord2_Union__0__ {
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
    float get_float() const;
    void set_float(const float& v);
    std::map<std::string, float > get_map() const;
    void set_map(const std::map<std::string, float >& v);
    _bigrecord2_Union__0__();
};

struct _bigrecord2_Union__1__ {
private:
    size_t idx_;
    boost::any value_;
public:
    size_t idx() const { return idx_; }
    float get_float() const;
    void set_float(const float& v);
    boost::array<uint8_t, 16> get_md5() const;
    void set_md5(const boost::array<uint8_t, 16>& v);
    _bigrecord2_Union__1__();
};

struct RootRecord {
    double mylong;
    int32_t anotherint;
    std::vector<uint8_t> bytes;
    Nested nestedrecord;
    std::map<std::string, int64_t > mymap;
    std::vector<double > myarray;
    ExampleEnum myenum;
    _bigrecord2_Union__0__ myunion;
    std::vector<uint8_t> anotherunion;
    Nested anothernested;
    bool newbool;
    _bigrecord2_Union__1__ myfixed;
};

inline
float _bigrecord2_Union__0__::get_float() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<float >(value_);
}

inline
void _bigrecord2_Union__0__::set_float(const float& v) {
    idx_ = 1;
    value_ = v;
}

inline
std::map<std::string, float > _bigrecord2_Union__0__::get_map() const {
    if (idx_ != 2) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::map<std::string, float > >(value_);
}

inline
void _bigrecord2_Union__0__::set_map(const std::map<std::string, float >& v) {
    idx_ = 2;
    value_ = v;
}

inline
float _bigrecord2_Union__1__::get_float() const {
    if (idx_ != 0) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<float >(value_);
}

inline
void _bigrecord2_Union__1__::set_float(const float& v) {
    idx_ = 0;
    value_ = v;
}

inline
boost::array<uint8_t, 16> _bigrecord2_Union__1__::get_md5() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<boost::array<uint8_t, 16> >(value_);
}

inline
void _bigrecord2_Union__1__::set_md5(const boost::array<uint8_t, 16>& v) {
    idx_ = 1;
    value_ = v;
}

inline _bigrecord2_Union__0__::_bigrecord2_Union__0__() : idx_(0) { }
inline _bigrecord2_Union__1__::_bigrecord2_Union__1__() : idx_(0), value_(float()) { }
}
namespace avro {
template<> struct codec_traits<testgen2::Nested> {
    static void encode(Encoder& e, const testgen2::Nested& v) {
        avro::encode(e, v.inval3);
        avro::encode(e, v.inval2);
        avro::encode(e, v.inval1);
    }
    static void decode(Decoder& d, testgen2::Nested& v) {
        avro::decode(d, v.inval3);
        avro::decode(d, v.inval2);
        avro::decode(d, v.inval1);
    }
};

template<> struct codec_traits<testgen2::ExampleEnum> {
    static void encode(Encoder& e, testgen2::ExampleEnum v) {
        e.encodeEnum(v);
    }
    static void decode(Decoder& d, testgen2::ExampleEnum& v) {
        v = static_cast<testgen2::ExampleEnum>(d.decodeEnum());
    }
};

template<> struct codec_traits<testgen2::_bigrecord2_Union__0__> {
    static void encode(Encoder& e, testgen2::_bigrecord2_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_float());
            break;
        case 2:
            avro::encode(e, v.get_map());
            break;
        }
    }
    static void decode(Decoder& d, testgen2::_bigrecord2_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 3) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                float vv;
                avro::decode(d, vv);
                v.set_float(vv);
            }
            break;
        case 2:
            {
                std::map<std::string, float > vv;
                avro::decode(d, vv);
                v.set_map(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen2::_bigrecord2_Union__1__> {
    static void encode(Encoder& e, testgen2::_bigrecord2_Union__1__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            avro::encode(e, v.get_float());
            break;
        case 1:
            avro::encode(e, v.get_md5());
            break;
        }
    }
    static void decode(Decoder& d, testgen2::_bigrecord2_Union__1__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            {
                float vv;
                avro::decode(d, vv);
                v.set_float(vv);
            }
            break;
        case 1:
            {
                boost::array<uint8_t, 16> vv;
                avro::decode(d, vv);
                v.set_md5(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen2::RootRecord> {
    static void encode(Encoder& e, const testgen2::RootRecord& v) {
        avro::encode(e, v.mylong);
        avro::encode(e, v.anotherint);
        avro::encode(e, v.bytes);
        avro::encode(e, v.nestedrecord);
        avro::encode(e, v.mymap);
        avro::encode(e, v.myarray);
        avro::encode(e, v.myenum);
        avro::encode(e, v.myunion);
        avro::encode(e, v.anotherunion);
        avro::encode(e, v.anothernested);
        avro::encode(e, v.newbool);
        avro::encode(e, v.myfixed);
    }
    static void decode(Decoder& d, testgen2::RootRecord& v) {
        avro::decode(d, v.mylong);
        avro::decode(d, v.anotherint);
        avro::decode(d, v.bytes);
        avro::decode(d, v.nestedrecord);
        avro::decode(d, v.mymap);
        avro::decode(d, v.myarray);
        avro::decode(d, v.myenum);
        avro::decode(d, v.myunion);
        avro::decode(d, v.anotherunion);
        avro::decode(d, v.anothernested);
        avro::decode(d, v.newbool);
        avro::decode(d, v.myfixed);
    }
};

}
#endif
