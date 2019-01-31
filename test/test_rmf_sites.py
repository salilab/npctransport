from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import RMF
import math
import IMP


class ConeTests(IMP.test.TestCase):
    def test_rmf_sites(self):
     print("p")
     RMF.set_log_level("Trace")
     seed = 1.0 # use time instead?
     IMP.random_number_generator.seed( seed )
      # create sim data and write to file then check
     config= IMP.npctransport.Configuration()
     IMP.npctransport.set_default_configuration(config)
     fg= IMP.npctransport.add_fg_type(config, type_name="fg0",
                                      number_of_beads=1,
                                      number=1,
                                      radius=10,
                                      interactions=2)
     config_name=self.get_tmp_file_name("rmf_sites.pb")
     cfg = open(config_name, "wb").write(config.SerializeToString())
     assignments_name= self.get_tmp_file_name("rmf_sites_assignments.pb")
     IMP.npctransport.assign_ranges(config_name, assignments_name,
                                    0, False, seed)
     rmf_name=self.get_tmp_file_name("rmf_sites.rmf")
     sd= IMP.npctransport.SimulationData(assignments_name, True)
     sd.set_rmf_file(rmf_name, False)
     wr= sd.get_rmf_sos_writer()
     wr.update_always()
     del wr
     del sd
     back= RMF.open_rmf_file_read_only(rmf_name)
     fg= back.get_root_node().get_children()[0].get_children()[0].get_children()[0]
     sites= fg.get_children()
     self.assertEqual(len(sites), 2)
     del back
     del cfg
     IMP.set_log_level(IMP.MEMORY)
     #del sd
     del fg
     for n in dir():
       del n
     print("done")


if __name__ == '__main__':
    IMP.test.main()
