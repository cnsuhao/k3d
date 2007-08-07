#ifndef TYPES_H
#define TYPES_h

#include <vector>

using std::vector;

namespace fluid_sim
{

class array3d_f : public vector<float>
{
public:
	array3d_f(int xcomps, int ycomps, int zcomps) :
		vector<float>(xcomps*ycomps*zcomps),
		m_xcomps(xcomps), m_ycomps(ycomps), m_zcomps(zcomps)
	{

	}

	float operator() (int i, int j, int k) const {
		return (*this)[k*m_xcomps*m_zcomps + m_xcomps*i + j];
	}
	


	float& operator() (int i, int j, int k) {
		return (*this)[k*m_xcomps*m_zcomps + m_xcomps*i + j];
	}


private:
	int m_xcomps, m_ycomps, m_zcomps;


};
	

}

#endif
