#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
using namespace boost;

class MarkovStatesMap{
  typedef std::set<IMP::ParticleType> ParticleTypeSet;
  typedef std::map<IMP::ParticleIndex, ParticleTypeSet> t_pi2ptype;

  struct Vertex{
    IMP::ParticleType particle_type;
  };
  typedef property<edge_weight_t, int> IntWeight;
  typedef property<edge_weight_t, double> DoubleWeight;
  typedef adjacency_list< vecS, vecS, directedS, IMP::ParticleType, IntWeight >
    DirectedGraph;

 private:
  mutable t_pi2ptype old_contacts;
  pi2ptypes contacts;

 public:

  void update_contact(IMP::Model m, IMP::ParticleIndexPair pip){
    IMP::Typed pt0(m, pip[0]) pt0;
    IMP::Typed pt1(m, pip[1]) pt1;

  }

  void update_contacts(
};
