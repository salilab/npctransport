/**
 *  \file csv.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/csv.h>
#include <IMP/algebra/GridD.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#ifndef IMP_NO_ACCUMULATORS
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#endif
#include <boost/algorithm/string/trim.hpp>
#include <boost/scoped_array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
namespace {
Strings get_split(char* buf) {
  Strings ret;
  std::istringstream iss(buf);
  std::string item;
  while (std::getline(iss, item, ',')) {
    boost::trim_if(item, boost::is_any_of("\" "));
    ret.push_back(item);
  }
  return ret;
}
}


CSV::CSV(base::TextInput in, const Strings &variables): Object(in.get_name()) {
  static const int bufsize=10000;
  boost::scoped_array<char> buf(new char[bufsize]);
  in.get_stream().getline(buf.get(), bufsize);
  IMP_USAGE_CHECK(buf[0]=='#', "No hash on first line");
  names_= get_split(buf.get());
  values_.resize(names_.size());
  do {
    in.get_stream().getline(buf.get(), bufsize);
    if (!in) break;
    Strings split= get_split(buf.get());
    if (split.size()!= names_.size()) {
      IMP_WARN("Can't parse line: " << buf.get() << " -> "
               << IMP::Showable(split) << " count " << names_.size()
               << std::endl);
      continue;
    }
    for (unsigned int i=0; i<split.size(); ++i) {
      double b= std::atof(split[i].c_str());
      values_[i].push_back(b);
    }
  } while(true);
  for (unsigned int i=0; i< variables.size(); ++i) {
    variables_.push_back(get_index(variables[i]));
  }
}


const Floats & CSV::get_grid_indexes(std::string name) const {
  for (unsigned int i=0; i< variables_.size(); ++i) {
    if (names_[variables_[i]]==name) {
      std::ostringstream oss;
      oss << "grid_index_"<< i;
      return values_[get_index(oss.str())];
    }
  }
  static Floats junk;
  return junk;
}
namespace {
typedef algebra::ExtendedGridIndexD<-1> Index;
Index get_hash_index(const Floats &v) {
  Ints seq(v.size());
  for (unsigned int i=0; i < seq.size(); ++i) {
    seq[i]= static_cast<int>(v[i]*1000);
  }
  return Index(seq.begin(), seq.end());
}
}

FloatsList CSV::get_maximum(const Strings& inputs,
                            const std::string &output,
                            const Strings &filters,
                            const Floats &filter_mins,
                            const Floats &filter_maxs) const {
  typedef compatibility::map<Index, double > Map;
  Map maxs;
  Ints indexes;
  Ints filter_indexes;
  int outindex= get_index(output);
  for (unsigned int i=0; i< inputs.size(); ++i) {
    indexes.push_back(get_index(inputs[i]));
  }
  for (unsigned int i=0; i< filters.size(); ++i) {
    filter_indexes.push_back(get_index(filters[i]));
  }
  for (unsigned int i=0; i< get_number_of_rows(); ++i) {
    Floats v(inputs.size());
    bool ok=true;
    for (unsigned int j=0; j< filters.size(); ++j) {
      if (values_[filter_indexes[j]][i] <= filter_maxs[j]
          && values_[filter_indexes[j]][i] >= filter_mins[j]) {

      } else {
        ok=false;
        break;
      }
    }
    if (!ok) continue;
    for (unsigned int j=0; j< inputs.size(); ++j) {
      v[j]= values_[indexes[j]][i];
    }
    Index vv=get_hash_index(v);

    if (maxs.find(vv)== maxs.end()) {
      maxs[vv]= values_[outindex][i];
    } else {
      maxs[vv] = std::max(maxs[vv], values_[outindex][i]);
    }
  }
  FloatsList ret; ret.reserve(maxs.size());
  for (Map::const_iterator
           it = maxs.begin(); it != maxs.end(); ++it){
    Floats v(inputs.size()+1);
    for (unsigned int i=0; i< inputs.size(); ++i) {
      v[i]= static_cast<double>(it->first[i])/1000.0;
    }
    v[v.size()-1]=it->second;
    ret.push_back(v);
  }
  return ret;
}


namespace {
void show_statistics( CSV *csv,
                     std::string field_name,
                      base::TextOutput out) {
#ifndef IMP_NO_ACCUMULATORS
  const Floats& data= csv->get_values(field_name);
  using namespace boost::accumulators;
  typedef features<tag::min,
      tag::max,
      tag::mean,
      tag::moment<2> > Features;
  typedef accumulator_set<double,Features > Accum;
  Accum acc;
  std::for_each(data.begin(), data.end(),
                boost::bind<void>(boost::ref(acc), _1));
  out.get_stream() << field_name << ", ";
  out.get_stream() << extract::min(acc) << ", ";
  out.get_stream() << extract::max(acc) << ", ";
  out.get_stream() << extract::mean(acc) << ", ";
  out.get_stream() << boost::accumulators::moment<2>(acc) << std::endl;
#endif
}
}

void show_statistics(CSV *csv,
                     const Strings& field_names,
                     base::TextOutput out) {
#ifndef IMP_NO_ACCUMULATORS
  base::OwnerPointer<CSV> ocsv(csv);
  out.get_stream() << "name, min, max, mean, stddev" << std::endl;
  for (unsigned int i=0; i< field_names.size(); ++i) {
    show_statistics(csv, field_names[i], out);
  }
#endif
}

namespace {
void show_joint_statistics( CSV *csv,
                     std::string field_name,
                            std::string output_name,
                      base::TextOutput out) {
#ifndef IMP_NO_ACCUMULATORS
  const Floats& data= csv->get_values(field_name);
  const Floats& outputdata= csv->get_values(output_name);
  using namespace boost::accumulators;
  typedef features<tag::covariance<double, tag::covariate1> > Features;
  typedef accumulator_set<double,Features > Accum;
  Accum acc;
  for (unsigned int i =0; i< data.size(); ++i) {
    acc(data[i], covariate1=outputdata[i]);
  }
  out.get_stream() << field_name << ", ";
  out.get_stream() << covariance(acc) << std::endl;
#endif
}
}

void show_joint_statistics(CSV *csv,
                     const Strings& field_names,
                           std::string output_name,
                     base::TextOutput out) {
#ifndef IMP_NO_ACCUMULATORS
  base::OwnerPointer<CSV> ocsv(csv);
  out.get_stream() << "name, covariance" << std::endl;
  for (unsigned int i=0; i< field_names.size(); ++i) {
    show_joint_statistics(csv, field_names[i], output_name, out);
  }
#endif
}


algebra::PrincipalComponentAnalysisKD
get_principal_components(CSV *csv,
                         const Strings& field_names) {
  base::OwnerPointer<CSV> ocsv(csv);
  algebra::VectorKDs vects(csv->get_number_of_rows());
  for (unsigned int i=0; i< vects.size(); ++i) {
    algebra::VectorKD kd= algebra::get_zero_vector_kd(field_names.size());
    for (unsigned int j=0; j< field_names.size(); ++j) {
      kd[j]=csv->get_values(field_names[j])[i];
    }
    vects[i]=kd;
  }
  return algebra::get_principal_components(vects);
}




void show_with_color(CSV *csv,
                     const Strings& base,
                     double ball_size,
                     std::string variable,
                     double min, double max,
                     Operation,
                     display::Writer *w) {
  base::OwnerPointer<CSV> ocsv(csv);
  base::OwnerPointer<display::Writer> ow(w);
  base::Vector<const Floats*> locations_indexes;
  for (unsigned int i=0; i< base.size(); ++i) {
    //location_values.push_back(&csv->get_values(locations[i]));
    locations_indexes.push_back(&csv->get_grid_indexes(base[i]));
  }
  const Floats &vs= csv->get_values(variable);
  typedef compatibility::map<algebra::ExtendedGridIndexKD, double> Map;
  Map values;
  for (unsigned int i=0; i< vs.size(); ++i) {
    Ints center;
    for (unsigned int j=0; j< base.size(); ++j) {
      center.push_back(locations_indexes[j]->at(i));
    }
    algebra::ExtendedGridIndexKD ind(center);
    if (values.find(ind) != values.end()) {
      values[ind]= std::max(values[ind], vs[i]);
    } else {
      values[ind]=vs[i];
    }
  }
  for (Map::const_iterator it= values.begin(); it != values.end(); ++it) {
    double svs= (it->second-min)/(max-min);
    if (svs>1) svs=1;
    if (svs<0) svs=0;
    if (svs >0) {
      IMP_NEW(display::SphereGeometry, sg,
              (algebra::Sphere3D(algebra::Vector3D(it->first.begin(),
                                                   it->first.end()),
                                 ball_size)));
      sg->set_color(display::get_jet_color(svs));
      w->add_geometry(sg);
    }
  }
}


IMPNPCTRANSPORT_END_NAMESPACE
