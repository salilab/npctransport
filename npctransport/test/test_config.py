import IMP
import IMP.test
import IMP.npctransport
import math

radius=5

class ConeTests(IMP.test.TestCase):
    def test_1(self):
        """Check discovery of ranges from config"""
        IMP.npctransport.show_ranges(self.get_input_file_name("test.config"))
    def test_2(self):
        """Check assign ranges"""
        IMP.npctransport.assign_ranges(self.get_input_file_name("test.config"),
                                       self.get_tmp_file_name("out.assign"), 1,
                                       True)
if __name__ == '__main__':
    IMP.test.main()
