import IMP
import IMP.test
import IMP.npctransport
import math

radius=5

class ConeTests(IMP.test.TestCase):
    def _make_config(self, name):
        config= IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        IMP.npctransport.create_range(config.interaction_k, .2, 20, 5)
        IMP.npctransport.create_range(config.backbone_k, .2, 20, 5)
        IMP.npctransport.create_range(config.nonspecific_range, 0, 5, 3, base=1)
        IMP.npctransport.create_range(config.nonspecific_k, .2, 20, 3)

        fg= IMP.npctransport.add_fg_type(config,
                                         number_of_beads=10,
                                         number=1,
                                         radius=10,
                                         interactions=1)
        IMP.npctransport.create_range(fg.rest_length_factor, .5, .8, 3)

        kap= IMP.npctransport.add_float_type(config,
                                             number=1,
                                             radius=10)
        IMP.npctransport.create_range(kap.number,0, 10,3)
        IMP.npctransport.create_range(kap.interactions, 1, 10, 3)

        interaction= IMP.npctransport.add_interaction(config, "fg0", "fg0")
        IMP.npctransport.create_range(interaction.is_on,0,1,2)

        interaction= IMP.npctransport.add_interaction(config, "fg0", "kap")
        f=open(name, "wb")
        f.write(config.SerializeToString())
    def test_1(self):
        """Check creating a configuration and assigning values"""

        config_name= self.get_tmp_file_name("configuration.pb")
        self._make_config(config_name)

        assignment_name= self.get_tmp_file_name("assignment.pb")
        IMP.npctransport.assign_ranges(config_name, assignment_name, 0, True)
        output= IMP.npctransport.Output()
        f=open(assignment_name, "rb")
        output.ParseFromString(f.read())
        assign= output.assignment
        self.assertAlmostEqual(assign.interaction_k.value, .2, delta=.0000001)
        self.assertEqual(assign.interactions[0].is_on.value, 0)
        num= IMP.npctransport.get_number_of_work_units(config_name)
        self.assertEqual(num, 5*5*3*3*3*3*3*2)
        assignment_name= self.get_tmp_file_name("final.pb")
        IMP.npctransport.assign_ranges(config_name, assignment_name, num-1, True)
        f=open(assignment_name, "rb")
        output.ParseFromString(f.read())
        assign= output.assignment
        self.assertEqual(assign.interaction_k.value, 20)
        self.assertEqual(assign.interactions[0].is_on.value, 1)
if __name__ == '__main__':
    IMP.test.main()
