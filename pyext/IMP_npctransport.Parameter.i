%include "IMP/npctransport/Parameter.h"
namespace IMP {
  namespace npctransport {
    %template(_DoubleParameter) Parameter< double >;
    %template(_IntParameter) Parameter< int >;
    %template(_BoolParameter) Parameter< bool >;
  }
}
