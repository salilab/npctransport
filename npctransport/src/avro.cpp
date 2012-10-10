#include <IMP/npctransport/avro.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

avro::ValidSchema get_avro_data_file_schema() {
  std::string path= get_data_path("avrodatafile.json");
  TextInput in(path);
  std::stringstream buf;
  buf << in.get_stream.rdbuf();
  avro::ValidSchema result= avro::compileJsonSchemaFromString(buf.str());
  return result;
}

IMPNPCTRANSPORT_END_NAMESPACE
