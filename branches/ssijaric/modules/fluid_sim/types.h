#ifndef TYPES_H
#define TYPES_h

#include <vector>

using std::vector;

namespace fluid_sim
{
template<class T>
class array3d : public vector<T>
{
public:
	array3d(int xcomps, int ycomps, int zcomps) :
		vector<T>(xcomps*ycomps*zcomps),
		m_xcomps(xcomps), m_ycomps(ycomps), m_zcomps(zcomps)
	{

	}

	T operator() (int i, int j, int k) const {
		return (*this)[k*m_xcomps*m_zcomps + m_xcomps*i + j];
	}
	


	T& operator() (int i, int j, int k) {
		return (*this)[k*m_xcomps*m_zcomps + m_xcomps*i + j];
	}

	int xcomps() { return m_xcomps; }
	int ycomps() { return m_ycomps; }
	int zcomps() { return m_zcomps; }


private:
	int m_xcomps, m_ycomps, m_zcomps;


};
	

}

#endif
