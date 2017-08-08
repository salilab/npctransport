from __future__ import print_function

import IMP
import IMP.test
import IMP.npctransport
import IMP
import math
from test_util import *


class InteractionsTests(IMP.test.TestCase):

    def _make_config_nonspecific_interactions(self,
                                              name,
                                              nonspecific_k= None,
                                              nonspecific_range= None,
                                              excluded_volume_k= None):
        config= IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        config.box_is_on.lower= 0
        config.nonspecific_k.lower= self.default_nonspecific_k
        config.nonspecific_range.lower= self.default_nonspecific_range
        config.excluded_volume_k.lower= self.default_excluded_volume_k
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
        if nonspecific_k is not None:
            i_interaction.nonspecific_k.lower= nonspecific_k
        if nonspecific_range is not None:
            i_interaction.nonspecific_range.lower= nonspecific_range
        if excluded_volume_k is not None:
            i_interaction.excluded_volume_k.lower= excluded_volume_k
        f=open(name, "wb")
        f.write(config.SerializeToString())

    def _test_score(self,
                    coords1,
                    coords2,
                    expected_score):
        #        self._test_score([0,0,0],[0,0,19],
        for bead in self.sd.get_beads():
            assert(IMP.core.XYZ.get_is_setup(bead))
            xyz= IMP.core.XYZ(bead)
            if IMP.core.Typed.get_is_setup(bead):
                t= IMP.core.Typed(bead)
                if t.get_type()==IMP.core.ParticleType("fg0"):
                    xyz.set_coordinates(coords1)
                if t.get_type()==IMP.core.ParticleType("inert0"):
                    xyz.set_coordinates(coords2)
#            print("Particle",bead,"xyz",xyz,"type",t)
        score= self.sd.get_bd().get_scoring_function().evaluate(False)
        self.assertAlmostEqual(score, expected_score)

    def test_interactions_from_protobuf(self):
        """
        Check that interactions are evaluated correctly
        when defined through protobuf
        """
        test_protobuf_installed(self)
        self.default_nonspecific_k= 100.0
        self.default_nonspecific_range=5.0
        self.default_excluded_volume_k= 5000.0
        seed = 1.0
        IMP.random_number_generator.seed(seed)
        config_name= self.get_tmp_file_name("configuration.pb")
#        config_name="config.pb"

        print("I. Testing with non-default values for non-specific interactions")
        nonspecific_k=0.1
        nonspecific_range= 10.0
        excluded_volume_k=20.0
        self._make_config_nonspecific_interactions(config_name,
                                                   nonspecific_k,
                                                   nonspecific_range,
                                                   excluded_volume_k)
        assignment_name= self.get_tmp_file_name("assignment.pb")
#        assignment_name="output.pb"
        print("Config file", config_name)
        print("Assignment file", assignment_name)
        IMP.npctransport.assign_ranges(config_name,
                                       assignment_name,
                                       0,
                                       True,
                                       seed)
        self.sd= IMP.npctransport.SimulationData(assignment_name,
                                                 False,
                                                 "")
        coords1=[0,0,0]
        for i in range(100):
            coords2= [0,0,i]
            expected_excluded_score= excluded_volume_k*(20.0-i) if i <= 20 else 0.0
            expected_attractive= -nonspecific_k*nonspecific_range
            if i>=20:
                expected_attractive= expected_attractive+nonspecific_k*(i-20.0)
            if expected_attractive>0.0:
                expected_attractive= 0.0
            expected_score= expected_excluded_score+expected_attractive
            self._test_score(coords1,
                             coords2,
                             expected_score)

        print("II. Testing with default values for non-specific interactions")
        self._make_config_nonspecific_interactions(config_name)
        assignment_name= self.get_tmp_file_name("assignment.pb")
#        assignment_name="output.pb"
        print("Config file", config_name)
        print("Assignment file", assignment_name)
        IMP.npctransport.assign_ranges(config_name,
                                       assignment_name,
                                       0,
                                       True,
                                       seed)
        self.sd= IMP.npctransport.SimulationData(assignment_name,
                                                 False,
                                                 "")
        coords1=[0,0,0]
        for i in range(100):
            coords2= [0,0,i]
            expected_excluded_score= self.default_excluded_volume_k*(20.0-i) if i <= 20 else 0.0
            expected_attractive= -self.default_nonspecific_k * self.default_nonspecific_range
            if i>=20:
                expected_attractive= expected_attractive + self.default_nonspecific_k*(i-20.0)
            if expected_attractive>0.0:
                expected_attractive= 0.0
            expected_score= expected_excluded_score+expected_attractive
            self._test_score(coords1,
                             coords2,
                             expected_score)





if __name__ == '__main__':
    IMP.test.main()
