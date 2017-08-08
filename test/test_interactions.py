from __future__ import print_function


import IMP
import IMP.test
import IMP.npctransport
import IMP
import math
from test_util import *

class InteractionsTests(IMP.test.TestCase):
    def _make_config_nonspecific_interactions(self, name):
        config= IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        IMP.npctransport.create_range(config.interaction_k, .2, 20, 5)
        IMP.npctransport.create_range(config.backbone_k, .2, 20, 5)
        IMP.npctransport.create_range(config.nonspecific_range, 0, 5, 3, base=1)
        IMP.npctransport.create_range(config.nonspecific_k, .2, 20, 3)

        fg= IMP.npctransport.add_fg_type(config,
                                         type_name="fg0",
                                         number_of_beads=1,
                                         number=1,
                                         radius=10,
                                         interactions=1)
        print("FG TYPE:", fg.type)
        fg.rest_length_factor.lower= 1.0
        inert= IMP.npctransport.add_float_type(config,
                                             number=1,
                                             radius=10,
                                             interactions=0,
                                             type_name="inert0")
        f_interaction= IMP.npctransport.add_interaction(config, "fg0", "fg0")
        i_interaction= IMP.npctransport.add_interaction(config, "fg0", "inert0",
                                                        interaction_k=1.0,
                                                        interaction_range=1.0)
        i_interaction.nonspecific_k.lower= 0.5
        i_interaction.nonspecific_range.lower= 10
        i_interaction.excluded_volume_k.lower= 5
        f=open(name, "wb")
        f.write(config.SerializeToString())
    def test_interactions_from_protobuf(self):
        """
        Check that interactions are evaluated correctly
        when defined through protobuf
        """
        test_protobuf_installed(self)
        seed = 1.0 # use time instead?
        IMP.random_number_generator.seed(seed)
        config_name= self.get_tmp_file_name("configuration.pb")
        self._make_config_nonspecific_interactions(config_name)
        assignment_name= self.get_tmp_file_name("assignment.pb")
        IMP.npctransport.assign_ranges(config_name, assignment_name, 0, True, seed)
        sd= IMP.npctransport.SimulationData(assignment_name, False, "")
        bead=sd.get_beads()
        for bead in sd.get_beads():
            assert(IMP.atom.XYZ.get_is_instance(bead))
            xyz= IMP.atom.XYZ(bead)
            if IMP.core.Typed.get_is_instance(bead):
                t= IMP.core.Typed(bead)
                if t.get_type()=="fg0":
                    xyz.set_coordinates([0,0,0])
                if t.get_type()=="inert0":
                    xyz.set_coordinates([0,0,25])
        score= sd.get_bd().get_scoring_function().evaluate(False)
        print("Score is {}".format(score))



if __name__ == '__main__':
    IMP.test.main()
