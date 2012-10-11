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


#ifndef TWEET_HH_4104945341__H_
#define TWEET_HH_4104945341__H_


#include "boost/any.hpp"
#include "Specific.hh"
#include "Encoder.hh"
#include "Decoder.hh"

namespace testgen3 {
struct _tweet_Union__0__ {
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
    int64_t get_long() const;
    void set_long(const int64_t& v);
    _tweet_Union__0__();
};

struct AvroPoint {
    double latitude;
    double longitude;
};

struct _tweet_Union__1__ {
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
    AvroPoint get_AvroPoint() const;
    void set_AvroPoint(const AvroPoint& v);
    _tweet_Union__1__();
};

struct _tweet_Union__2__ {
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
    std::string get_string() const;
    void set_string(const std::string& v);
    _tweet_Union__2__();
};

struct AvroDateTime {
    std::string dateTimeString;
};

struct _tweet_Union__3__ {
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
    std::string get_string() const;
    void set_string(const std::string& v);
    _tweet_Union__3__();
};

struct AvroKnowableOptionString {
    bool known;
    _tweet_Union__3__ data;
};

struct AvroKnowableListString {
    bool known;
    std::vector<std::string > data;
};

struct AvroKnowableBoolean {
    bool known;
    bool data;
};

struct _tweet_Union__4__ {
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
    AvroPoint get_AvroPoint() const;
    void set_AvroPoint(const AvroPoint& v);
    _tweet_Union__4__();
};

struct AvroKnowableOptionPoint {
    bool known;
    _tweet_Union__4__ data;
};

struct AvroTweetMetadata {
    AvroKnowableOptionString inReplyToScreenName;
    AvroKnowableListString mentionedScreenNames;
    AvroKnowableListString links;
    AvroKnowableListString hashtags;
    AvroKnowableBoolean isBareCheckin;
    AvroKnowableBoolean isBareRetweet;
    AvroKnowableBoolean isRetweet;
    AvroKnowableOptionString venueID;
    AvroKnowableOptionPoint venuePoint;
};

struct AvroTweet {
    int64_t ID;
    std::string text;
    std::string authorScreenName;
    std::string authorProfileImageURL;
    _tweet_Union__0__ authorUserID;
    _tweet_Union__1__ location;
    _tweet_Union__2__ placeID;
    AvroDateTime createdAt;
    AvroTweetMetadata metadata;
};

inline
int64_t _tweet_Union__0__::get_long() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<int64_t >(value_);
}

inline
void _tweet_Union__0__::set_long(const int64_t& v) {
    idx_ = 1;
    value_ = v;
}

inline
AvroPoint _tweet_Union__1__::get_AvroPoint() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<AvroPoint >(value_);
}

inline
void _tweet_Union__1__::set_AvroPoint(const AvroPoint& v) {
    idx_ = 1;
    value_ = v;
}

inline
std::string _tweet_Union__2__::get_string() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::string >(value_);
}

inline
void _tweet_Union__2__::set_string(const std::string& v) {
    idx_ = 1;
    value_ = v;
}

inline
std::string _tweet_Union__3__::get_string() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<std::string >(value_);
}

inline
void _tweet_Union__3__::set_string(const std::string& v) {
    idx_ = 1;
    value_ = v;
}

inline
AvroPoint _tweet_Union__4__::get_AvroPoint() const {
    if (idx_ != 1) {
        throw avro::Exception("Invalid type for union");
    }
    return boost::any_cast<AvroPoint >(value_);
}

inline
void _tweet_Union__4__::set_AvroPoint(const AvroPoint& v) {
    idx_ = 1;
    value_ = v;
}

inline _tweet_Union__0__::_tweet_Union__0__() : idx_(0) { }
inline _tweet_Union__1__::_tweet_Union__1__() : idx_(0) { }
inline _tweet_Union__2__::_tweet_Union__2__() : idx_(0) { }
inline _tweet_Union__3__::_tweet_Union__3__() : idx_(0) { }
inline _tweet_Union__4__::_tweet_Union__4__() : idx_(0) { }
}
namespace avro {
template<> struct codec_traits<testgen3::_tweet_Union__0__> {
    static void encode(Encoder& e, testgen3::_tweet_Union__0__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_long());
            break;
        }
    }
    static void decode(Decoder& d, testgen3::_tweet_Union__0__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                int64_t vv;
                avro::decode(d, vv);
                v.set_long(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen3::AvroPoint> {
    static void encode(Encoder& e, const testgen3::AvroPoint& v) {
        avro::encode(e, v.latitude);
        avro::encode(e, v.longitude);
    }
    static void decode(Decoder& d, testgen3::AvroPoint& v) {
        avro::decode(d, v.latitude);
        avro::decode(d, v.longitude);
    }
};

template<> struct codec_traits<testgen3::_tweet_Union__1__> {
    static void encode(Encoder& e, testgen3::_tweet_Union__1__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_AvroPoint());
            break;
        }
    }
    static void decode(Decoder& d, testgen3::_tweet_Union__1__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                testgen3::AvroPoint vv;
                avro::decode(d, vv);
                v.set_AvroPoint(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen3::_tweet_Union__2__> {
    static void encode(Encoder& e, testgen3::_tweet_Union__2__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_string());
            break;
        }
    }
    static void decode(Decoder& d, testgen3::_tweet_Union__2__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                std::string vv;
                avro::decode(d, vv);
                v.set_string(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen3::AvroDateTime> {
    static void encode(Encoder& e, const testgen3::AvroDateTime& v) {
        avro::encode(e, v.dateTimeString);
    }
    static void decode(Decoder& d, testgen3::AvroDateTime& v) {
        avro::decode(d, v.dateTimeString);
    }
};

template<> struct codec_traits<testgen3::_tweet_Union__3__> {
    static void encode(Encoder& e, testgen3::_tweet_Union__3__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_string());
            break;
        }
    }
    static void decode(Decoder& d, testgen3::_tweet_Union__3__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                std::string vv;
                avro::decode(d, vv);
                v.set_string(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen3::AvroKnowableOptionString> {
    static void encode(Encoder& e, const testgen3::AvroKnowableOptionString& v) {
        avro::encode(e, v.known);
        avro::encode(e, v.data);
    }
    static void decode(Decoder& d, testgen3::AvroKnowableOptionString& v) {
        avro::decode(d, v.known);
        avro::decode(d, v.data);
    }
};

template<> struct codec_traits<testgen3::AvroKnowableListString> {
    static void encode(Encoder& e, const testgen3::AvroKnowableListString& v) {
        avro::encode(e, v.known);
        avro::encode(e, v.data);
    }
    static void decode(Decoder& d, testgen3::AvroKnowableListString& v) {
        avro::decode(d, v.known);
        avro::decode(d, v.data);
    }
};

template<> struct codec_traits<testgen3::AvroKnowableBoolean> {
    static void encode(Encoder& e, const testgen3::AvroKnowableBoolean& v) {
        avro::encode(e, v.known);
        avro::encode(e, v.data);
    }
    static void decode(Decoder& d, testgen3::AvroKnowableBoolean& v) {
        avro::decode(d, v.known);
        avro::decode(d, v.data);
    }
};

template<> struct codec_traits<testgen3::_tweet_Union__4__> {
    static void encode(Encoder& e, testgen3::_tweet_Union__4__ v) {
        e.encodeUnionIndex(v.idx());
        switch (v.idx()) {
        case 0:
            e.encodeNull();
            break;
        case 1:
            avro::encode(e, v.get_AvroPoint());
            break;
        }
    }
    static void decode(Decoder& d, testgen3::_tweet_Union__4__& v) {
        size_t n = d.decodeUnionIndex();
        if (n >= 2) { throw avro::Exception("Union index too big"); }
        switch (n) {
        case 0:
            d.decodeNull();
            v.set_null();
            break;
        case 1:
            {
                testgen3::AvroPoint vv;
                avro::decode(d, vv);
                v.set_AvroPoint(vv);
            }
            break;
        }
    }
};

template<> struct codec_traits<testgen3::AvroKnowableOptionPoint> {
    static void encode(Encoder& e, const testgen3::AvroKnowableOptionPoint& v) {
        avro::encode(e, v.known);
        avro::encode(e, v.data);
    }
    static void decode(Decoder& d, testgen3::AvroKnowableOptionPoint& v) {
        avro::decode(d, v.known);
        avro::decode(d, v.data);
    }
};

template<> struct codec_traits<testgen3::AvroTweetMetadata> {
    static void encode(Encoder& e, const testgen3::AvroTweetMetadata& v) {
        avro::encode(e, v.inReplyToScreenName);
        avro::encode(e, v.mentionedScreenNames);
        avro::encode(e, v.links);
        avro::encode(e, v.hashtags);
        avro::encode(e, v.isBareCheckin);
        avro::encode(e, v.isBareRetweet);
        avro::encode(e, v.isRetweet);
        avro::encode(e, v.venueID);
        avro::encode(e, v.venuePoint);
    }
    static void decode(Decoder& d, testgen3::AvroTweetMetadata& v) {
        avro::decode(d, v.inReplyToScreenName);
        avro::decode(d, v.mentionedScreenNames);
        avro::decode(d, v.links);
        avro::decode(d, v.hashtags);
        avro::decode(d, v.isBareCheckin);
        avro::decode(d, v.isBareRetweet);
        avro::decode(d, v.isRetweet);
        avro::decode(d, v.venueID);
        avro::decode(d, v.venuePoint);
    }
};

template<> struct codec_traits<testgen3::AvroTweet> {
    static void encode(Encoder& e, const testgen3::AvroTweet& v) {
        avro::encode(e, v.ID);
        avro::encode(e, v.text);
        avro::encode(e, v.authorScreenName);
        avro::encode(e, v.authorProfileImageURL);
        avro::encode(e, v.authorUserID);
        avro::encode(e, v.location);
        avro::encode(e, v.placeID);
        avro::encode(e, v.createdAt);
        avro::encode(e, v.metadata);
    }
    static void decode(Decoder& d, testgen3::AvroTweet& v) {
        avro::decode(d, v.ID);
        avro::decode(d, v.text);
        avro::decode(d, v.authorScreenName);
        avro::decode(d, v.authorProfileImageURL);
        avro::decode(d, v.authorUserID);
        avro::decode(d, v.location);
        avro::decode(d, v.placeID);
        avro::decode(d, v.createdAt);
        avro::decode(d, v.metadata);
    }
};

}
#endif
