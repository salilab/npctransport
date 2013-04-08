import IMP
import IMP.test
import IMP.npctransport
import RMF
import math
import IMP.base


class ConeTests(IMP.test.TestCase):
    def test_rmf_sites(self):
     print "p"
     RMF.set_log_level("Trace")
     seed = 1.0 # use time instead?
     IMP.base.random_number_generator.seed( seed )
      # create sim data and write to file then check
     config= IMP.npctransport.Configuration()
     IMP.npctransport.set_default_configuration(config)
     fg= IMP.npctransport.add_fg_type(config,
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
     sd= IMP.npctransport.SimulationData(assignments_name, True, rmf_name)
     print "hi"
     wr= sd.get_rmf_sos_writer()
     print "ho"
     wr.update_always()
     del wr
     print "hoo"
     del sd
     back= RMF.open_rmf_file_read_only(rmf_name)
     print "hey"
     fg= back.get_root_node().get_children()[0].get_children()[0].get_children()[0]
     print "hee"
     sites= fg.get_children()
     print "heeya"
     self.assertEqual(len(sites), 2)
     del back
     del cfg
     print "hoohaa"
     IMP.base.set_log_level(IMP.base.MEMORY)
     #del sd
     del fg
     for n in dir():
       del n
     print "done"
if __name__ == '__main__':
  try:
    IMP.test.main()
  except:
    print "oops"
