/**
 *  \file util.h
 *  \brief Misceleneous internal utility functions
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_UTIL_H
#define IMPNPCTRANSPORT_INTERNAL_UTIL_H

#include "../npctransport_config.h"
#include <functional>
#include <iterator>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

//! from boost - just to support various problematic distributions
template<class Iterator, class Comp>
  inline Iterator is_sorted_until (Iterator first, Iterator last, Comp c) {
  if (first == last)
    return last;

  Iterator it = first; ++it;

  for (; it != last; first = it, ++it)
    if (c(*it, *first))
      return it;

  return it;
}

//! from boost - just to support various problematic distributions
template<class Iterator>
inline Iterator is_sorted_until (Iterator first, Iterator last) {
  typedef typename std::iterator_traits<Iterator>::value_type
    value_type;

  typedef std::less<value_type> c;

  return is_sorted_until(first, last, c());
}

//! from boost - just to support various problematic distributions
template<class Iterator, class Comp>
  inline bool is_sorted (Iterator first, Iterator last, Comp c) {
  return is_sorted_until(first, last, c) == last;
}

//! from boost - just to support various problematic distributions
template<class Iterator>
inline bool is_sorted (Iterator first, Iterator last) {
  return is_sorted_until(first, last) == last;
}




IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_UTIL_H */
