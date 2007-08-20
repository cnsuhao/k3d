#ifndef TYPES_H
#define TYPES_H

#include <vector>

using std::vector;

namespace fluid_sim
{

struct idx {
	idx(int a,int b,int c) : i(a), j(b), k(c) { }
	int i, j, k;
};

class particle {
public:
	particle(const k3d::point3& p) : m_pos(p) { }
private:
	k3d::point3 m_pos;
	k3d::point3 m_vel;
};


}

#endif
