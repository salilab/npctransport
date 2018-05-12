from __future__ import print_function


import IMP
import IMP.test
import IMP.npctransport
import IMP
import math
from test_util import *
radius=5

class ProtobufTests(IMP.test.TestCase):
    def _make_config(self,
                     name,
                     interaction_range= 1,
                     nonspecific_range= 1,
                     fg_radii= [10, 100],
                     kap_radii= [50, 500],
                     slack= 2.5):
        config= IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        config.interaction_k.lower = 1
        config.backbone_k.lower= 1
        config.nonspecific_k.lower= 1 #arbitrary
        config.nonspecific_range.lower= nonspecific_range
        config.slack.lower= slack
        fgs= []
        for fg_radius in fg_radii:
            fg= IMP.npctransport.add_fg_type \
                ( config,
                  type_name="fg{0}".format(fg_radius),
                  radius=fg_radius,
                  number_of_beads= 10,
                  number= 1,
                  interactions= 1 )
            fgs.append(fg)
        kaps= []
        for kap_radius in kap_radii:
            kap= IMP.npctransport.add_float_type \
                 ( config,
                   type_name="kap{0}".format(kap_radius),
                   radius=kap_radius,
                   number= 1,
                   interactions= 4 )
            kaps.append(kap)
        for a in (fgs+kaps):
            for b in (fgs+kaps):
                interaction= IMP.npctransport.add_interaction \
                             ( config,
                               a.type, b.type,
                               interaction_range= interaction_range,
                               interaction_k= 1 )
        f=open(name, "wb")
        f.write(config.SerializeToString())
    def _get_assignment(self, assignment_filename):
        seed = 1.0
        IMP.random_number_generator.seed(seed)
        config_filename= self.get_tmp_file_name("configuration.pb")
        self._make_config(config_filename)
        IMP.npctransport.assign_ranges(config_filename,
                                       assignment_filename, 0, True, seed)
        output= IMP.npctransport.Output()
        a= output.assignment
        f=open(assignment_filename, "rb")
        output.ParseFromString(f.read())
        return output.assignment
    def test_close_range_params(self):
        """Check creating a configuration and assigning values"""
        test_protobuf_installed(self)
        assignment_filename= self.get_tmp_file_name("assignment.pb")
        a= self._get_assignment(assignment_filename)
if __name__ == '__main__':
    IMP.test.main()
