import IMP
import IMP.test
import IMP.npctransport
import RMF
import math


class ConeTests(IMP.test.TestCase):
    def test_rmf_sites(self):
      # create sim data and write to file then check
      config= IMP.npctransport.Configuration()
      IMP.npctransport.set_default_configuration(config)
      fg= IMP.npctransport.add_fg_type(config,
                                       number_of_beads=1,
                                       number=1,
                                       radius=10,
                                       interactions=2)
      config_name=self.get_tmp_file_name("rmf_sites.pb")
      open(config_name, "wb").write(config.SerializeToString())
      assignments_name= self.get_tmp_file_name("rmf_sites_assignments.pb")
      IMP.npctransport.assign_ranges(config_name, assignments_name,
                                     0, False)
      sd= IMP.npctransport.SimulationData(assignments_name, True)
      rmf_name=self.get_tmp_file_name("rmf_sites.rmf")

      wr= sd.create_rmf_writer(rmf_name)
      wr.update_always()
      back= RMF.open_rmf_file_read_only(rmf_name)
      fg= back.get_root_node().get_children()[0].get_children()[0].get_children()[0]
      sites= fg.get_children()
      self.assertEqual(len(sites), 2)
if __name__ == '__main__':
    IMP.test.main()
