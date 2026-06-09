import IMP
import IMP.test
import IMP.npctransport
import IMP.display


def optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords, n_optimizations=100, verbose=False):
    # setup model
    m = IMP.Model()
    p = IMP.Particle(m)
    d = IMP.core.XYZR.setup_particle(p)
    d.set_coordinates_are_optimized(True)
    d.set_radius(particle_radius)
    bb = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                    IMP.algebra.Vector3D(boxw, boxw, boxw))
    zbss = IMP.npctransport.ZBiasSingletonScore(z_bias_k, z_bias_max_r, z_bias_z)
    sr = IMP.core.SingletonRestraint(m, zbss, p.get_index(), "zbias")
    d.set_coordinates(particle_init_coords)
    cg = IMP.core.SteepestDescent(m)
    cg.set_scoring_function(sr)
    cg.set_log_level(IMP.VERBOSE)
    for _ in range(n_optimizations):
        score = cg.optimize(3)
        if verbose:
            print(d.get_coordinates(), score)
    return d.get_x(), d.get_y(), d.get_z()


class ZBiasTests(IMP.test.TestCase):
    def test_z_bias_1(self):
        """Check z-axis bias singleton score - particle hovers around z=0 (from top)"""

        # setup parameters
        particle_radius=1
        z_bias_k = 10.0
        z_bias_z = 0
        z_bias_max_r = 1000.0
        boxw= 20
        particle_init_coords = [1, 1, 15]

        x, y, z = optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertTrue((z - z_bias_z) <= 0.1 and
                        (z - z_bias_z) >= -0.1 and
                        x == particle_init_coords[0] and
                        y == particle_init_coords[1])

    def test_z_bias_2(self):
        """Check z-axis bias singleton score - particle hovers around z=0 (from bottom)"""

        # setup parameters
        particle_radius=1
        z_bias_k = 10.0
        z_bias_z = 0
        z_bias_max_r = 1000.0
        boxw= 20
        particle_init_coords = [1, 1, -15]

        x, y, z = optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertTrue((z - z_bias_z) <= 0.1 and
                        (z - z_bias_z) >= -0.1 and
                        x == particle_init_coords[0] and
                        y == particle_init_coords[1])

    def test_z_bias_3(self):
        """Check z-axis bias singleton score - particle hovers around z=3 (from top)"""

        # setup parameters
        particle_radius=1
        z_bias_k = 10.0
        z_bias_z = 3
        z_bias_max_r = 1000.0
        boxw= 20
        particle_init_coords = [1, 1, 15]

        x, y, z = optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertTrue((z - z_bias_z) <= 0.1 and
                        (z - z_bias_z) >= -0.1 and
                        x == particle_init_coords[0] and
                        y == particle_init_coords[1])

    def test_z_bias_4(self):
        """Check z-axis bias singleton score - particle hovers around z=3 (from bottom)"""

        # setup parameters
        particle_radius=1
        z_bias_k = 10.0
        z_bias_z = 3
        z_bias_max_r = 1000.0
        boxw= 20
        particle_init_coords = [1, 1, -15]

        x, y, z = optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertTrue((z - z_bias_z) <= 0.1 and
                        (z - z_bias_z) >= -0.1 and
                        x == particle_init_coords[0] and
                        y == particle_init_coords[1])

    def test_z_bias_5(self):
        """Check z-axis bias singleton score - lower k slower convergence than higher k"""

        # setup parameters
        particle_radius=1
        z_bias_k_1 = 0.1
        z_bias_k_2 = 10.0
        z_bias_z = 0
        z_bias_max_r = 1000.0
        boxw= 20
        particle_init_coords = [1, 1, 15]

        x1, y1, z1 = optimize_particle_z_bias(particle_radius, z_bias_k_1, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        x2, y2, z2 = optimize_particle_z_bias(particle_radius, z_bias_k_2, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertFalse((z1 - z_bias_z) <= 0.1 and
                        (z1 - z_bias_z) >= -0.1)
        self.assertTrue((z2 - z_bias_z) <= 0.1 and
                        (z2 - z_bias_z) >= -0.1 and
                        x2 == particle_init_coords[0] and
                        y2 == particle_init_coords[1] and
                        z1 != particle_init_coords[2])

    def test_z_bias_6(self):
        """Check z-axis bias singleton score - particle outside max_r does not move"""

        # setup parameters
        particle_radius=1
        z_bias_k = 10.0
        z_bias_z = 0
        z_bias_max_r = 5
        boxw= 20
        particle_init_coords = [10, 10, 15]

        x, y, z = optimize_particle_z_bias(particle_radius, z_bias_k, z_bias_z, z_bias_max_r, boxw, particle_init_coords)
        self.assertTrue(x == particle_init_coords[0] and
                        y == particle_init_coords[1] and
                        z == particle_init_coords[2])


if __name__ == '__main__':
    IMP.test.main()
