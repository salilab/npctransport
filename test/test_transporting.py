from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math
from test_util import *

radius=20

class Tests(IMP.test.TestCase):
    def test_transporting(self):
        m=IMP.Model()
        p=create_diffusing_rb_particle(m,radius)
        t=IMP.npctransport.Transporting.setup_particle(p)
        t.set_is_last_entry_from_top(True);
        self.assertTrue(t.get_is_last_entry_from_top())
        t.set_is_last_entry_from_top(False);
        self.assertFalse(t.get_is_last_entry_from_top())
        t.set_last_tracked_z(50.0);
        self.assertEqual(t.get_last_tracked_z(),50.0)
        t.set_n_entries_bottom(10)
        self.assertEqual(t.get_n_entries_bottom(),10)
        t.set_n_entries_top(20)
        self.assertEqual(t.get_n_entries_top(),20)
        # test decoration
        t2=IMP.npctransport.Transporting(p)
        self.assertEqual(t.get_particle_index(), t2.get_particle_index())


if __name__ == '__main__':
    IMP.test.main()
