#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/avro.h>
#include <ValidSchema.hh>
#include <Compiler.hh>
#include <Stream.hh>
#include <IMP/base/file.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

avro::ValidSchema get_avro_data_file_schema() {
  std::string path = get_data_path("AvroDataFileData.json");
  std::auto_ptr<avro::InputStream> is = avro::fileInputStream(path.c_str());
  return avro::compileJsonSchemaFromStream(*is);
}

IMPNPCTRANSPORT_END_NAMESPACE
