Import('env')
import scons_tools.dependency.boost
import scons_tools.install
import scons_tools.paths
from SCons.Script import Copy

import os

if env["IMP_PASS"]=="BUILD":
  # what ridiculousness scons makes us do

  path=Dir("#/build/data/npctransport").abspath
  proto=path+"/npctransport.proto"
  #print proto, path
  # work around dumb script
  try:
    os.makedirs(Dir("#/build/lib/IMP/npctransport/").abspath)
  except:
    pass
  try:
    os.makedirs(Dir("#/build/src/npctransport").abspath)
  except:
    pass
  proto_files = env.Protoc(
      [],
      proto,
      PROTOCPROTOPATH=[path],
      PROTOCPYTHONOUTDIR=Dir("#/build/lib/IMP/npctransport/").abspath,
      PROTOCOUTDIR = Dir("#/build/src/npctransport/").abspath,
      PROTOCINCLUDE = "IMP/npctransport/npctransport_config.h",
      PROTOCCPPOUTFLAGS = "dllexport_decl=IMPNPCTRANSPORTEXPORT:"
      )
  env.Command(File("#/build/include/IMP/npctransport/internal/npctransport.pb.h"),
              File("#/build/src/npctransport/npctransport.pb.h"),
              Copy("$TARGET", "$SOURCE"))
  env.Command(File("#/build/src/npctransport/npctransport.pb.cpp"),
              File("#/build/src/npctransport/npctransport.pb.cc"),
              Copy("$TARGET", "$SOURCE"))

  def build_avro(target, source, env):
        cmd="avrogencpp -i %(json)s -o %(header)s -n IMP_npctransport"%{
     "json":source[0].abspath,
     "header":target[0].abspath}
        print "cmd=", cmd
        env.Execute(cmd)
  env['BUILDERS']['AvroCpp']=Builder(action=build_avro)
  env.AvroCpp(source=[File("data/AvroDataFileData.json")],
              target=[File("#/build/include/IMP/npctransport/AvroDataFileData.h")])


oenv=env.IMPModuleBuild()
