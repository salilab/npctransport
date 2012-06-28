/**
 *  \file csv.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CSV_H
#define IMPNPCTRANSPORT_CSV_H

#include "npctransport_config.h"
#include <IMP/base/file.h>
#include <IMP/algebra/eigen_analysis.h>
#include <IMP/base/Object.h>
#include <IMP/display/Writer.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
/** Take a set of rows of csv data structured in a grid. The data
    is stored based on that grid.
 */
class IMPNPCTRANSPORTEXPORT CSV: public base::Object {
  base::Vector<std::string> names_;
  typedef algebra::SparseGridStorageD<-1, algebra::VectorKD,
      algebra::UnboundedGridRangeKD> Storage;
  typedef algebra::GridD<-1, Storage, algebra::VectorKD> Grid;
  Grid grid;
  Ints variables_;
  base::Vector<Floats> values_;
  unsigned int get_index(std::string name) const {
    return std::find(names_.begin(),
                     names_.end(), name)-names_.begin();
  }
  // return the Index for the current row, identifying the parameter values
  algebra::ExtendedGridIndexKD get_row_index(unsigned int i) const;
 public:
  // average over rows with identical variables
  CSV(base::TextInput in, const Strings& variables);
  const Floats& get_values(std::string name) const {
    return values_[get_index(name)];
  }
  unsigned int get_number_of_rows() const {
    return values_[0].size();
  }
  const Floats& get_grid_indexes(std::string name) const;
  /** Return the maximum value for a given variable,
      over each set of values for the input variables.*/
  FloatsList get_maximum(const Strings& inputs,
                         const std::string &output,
                         const Strings &filters=Strings(),
                         const Floats &filter_mins=Floats(),
                         const Floats &filter_maxs=Floats()) const;
  IMP_OBJECT_INLINE(CSV,IMP_UNUSED(out),);
};

IMPNPCTRANSPORTEXPORT
void show_statistics( CSV *csv,
                      const Strings& field_names,
                      base::TextOutput out);

IMPNPCTRANSPORTEXPORT
void show_joint_statistics( CSV *csv,
                             const Strings& field_names,
                             std::string output_name,
                             base::TextOutput out);

IMPNPCTRANSPORTEXPORT
algebra::PrincipalComponentAnalysisKD
get_principal_components(CSV *csv,
                         const Strings& field_names);

enum Operation {MIN, MAX};
IMPNPCTRANSPORTEXPORT
void show_with_color(CSV *csv,
                     const Strings& base,
                     double ball_size,
                     std::string variable,
                     double min, double max,
                     Operation op,
                     display::Writer *w);


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_CSV_H */
