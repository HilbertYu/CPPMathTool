#ifndef MFUNCTION_H
#define MFUNCTION_H

#include "MVector.h"

template <size_t Dim, class scale_t>
class MFunction
{
    scale_t m_domain_h;
    mutable std::map< MVector<scale_t>, scale_t> m_func_val;
public:

    MFunction(void):
        m_domain_h(1),
        m_func_val()
    {}

    virtual scale_t impFunc(const MVector<scale_t> & x) const = 0;

    virtual scale_t func(const MVector<scale_t> & x) const
    {
        typedef std::map< MVector<scale_t>, scale_t> map_t;
        typename map_t::const_iterator itr = m_func_val.find(x);

        if (itr != m_func_val.end())
            return itr->second;

        return impFunc(x);
    }

    scale_t parital_dif(int dir, const MVector<scale_t> & x) const
    {
        MVector<scale_t> rx = x;
        MVector<scale_t> lx = x;
        rx[dir] += m_domain_h;
        lx[dir] -= m_domain_h;

        return (func(rx) - func(lx))/(2*m_domain_h);
    }

    MVector<scale_t> gradient(const MVector<scale_t> & x) const
    {
        MVector<scale_t> ret;
        for (size_t i = 0; i != Dim; ++i)
        {
            ret.push_back(parital_dif(i, x));
        }
        return ret;
    }

    virtual ~MFunction(void)
    {}

};

#endif /* end of include guard: MFUNCTION_H */
