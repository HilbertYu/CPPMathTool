#ifndef MVECTOR_H
#define MVECTOR_H

#include <stdio.h>
#include <vector>
#include <deque>
#include <iostream>
#include <assert.h>
#include <math.h>

template <typename T>
class MVector
{
    typedef std::deque<T> vector_t;

    vector_t m_vec;
public:

    enum OperatorType
    {
        OP_PLUS,
        OP_MINUS,
        OP_MULT,
        OP_DIVIDE,

    };

    MVector(void):
        m_vec()
    {

    }

    template <typename SrcT>
    MVector(SrcT itr,  SrcT end):
        m_vec()
    {
        attach(itr, end);
    }

    template <typename SrcT>
    MVector(SrcT x):
        m_vec()
    {
        attach(x.begin(), x.end());
    }

    void zeros(size_t n_z)
    {
        m_vec.clear();
        
        for (size_t i = 0; i < n_z; ++i)
            m_vec.push_back(0);
    }

    size_t size(void ) const
    {
        return m_vec.size();
    }

    bool empty(void) const
    {
        return m_vec.empty();
    }

    void clear(void) const
    {
        m_vec.clear();
    }

    void push_back(const T & v)
    {
        m_vec.push_back(v);
    }

    void push_front(const T & v)
    {
        m_vec.push_front(v);
    }

    T pop_back(void)
    {
        assert(!empty());
        const T & ret = m_vec.back();
        m_vec.pop_back();
        return ret;
    }

    T pop_front(const T & v)
    {
        assert(!empty());
        const T & ret = m_vec.back();
        m_vec.pop_front();
        return ret;
    }
    //=====================

    T sum(void) const
    {
        assert(!empty());
        T ret = 0;
        typename vector_t::const_iterator itr = m_vec.begin();
        for (; itr != m_vec.end(); ++itr)
            ret += *itr;
        return ret;
    }

    T mean(void) const
    {
        assert(!empty());
        return sum()/size();
    }

    T var(int nl = 0) const
    {
        assert(size() - nl > 0);
        T ret = 0;

        typename vector_t::const_iterator itr = m_vec.begin();
        T m = mean();
        for (; itr != m_vec.end(); ++itr)
            ret += pow(*itr - m, 2);
        return ret/(size()-nl);
    }

    T stddev(int nl = 0) const
    {
        assert(!empty());
        T v = var(nl);
        return sqrt(v);
    }

    T norm(void) const
    {
        return sqrt((*this**this).sum());
    }

    void eachPower(int n)
    {
        assert(n >= 0);

        typename vector_t::iterator itr = m_vec.begin();
        for (; itr != m_vec.end(); ++itr)
        {
            *itr = pow(*itr, n);
        }

    }

    //======operators===========

    T &operator[](int idx)
    {
        return m_vec.at(idx);
    }

    const T &operator[](int idx) const
    {
        return m_vec.at(idx);
    }


    //----- imp ---------

    template <class Op1_t, class Op2_t>
    static T impBinOperator(const Op1_t & op1, const Op2_t & op2, OperatorType op)
    {
        T ret;

        switch (op)
        {
            case OP_PLUS:
                ret = op1 + op2;
                break;
            case OP_MINUS:
                ret = op1 - op2;
                break;
            case OP_MULT:
                ret = op1 * op2;
                break;
            case OP_DIVIDE:
                ret = op1 / op2;
                break;
            default:
                assert(false);
        }
        return ret;
    }


    static
    MVector opVecVec(const MVector & v1, const MVector & v2, OperatorType op)
    {
        assert(v1.size() == v2.size());
        MVector ret;

        for (size_t i = 0; i < v1.size(); ++i)
        {
            ret.push_back(impBinOperator(v1[i], v2[i], op));
        }

        return ret;

    }

    template <class scale_t>
    static MVector opVecScale(const MVector & v, const scale_t & s, OperatorType op)
    {
        assert(!v.empty());
        MVector ret;

        for (size_t i = 0; i != v.size(); ++i)
        {
            ret.push_back(impBinOperator(v[i], s, op));
        }

        return ret;
    }


    //-------------------
    bool operator<(const MVector &v) const
    {
        return m_vec < v.m_vec;
    }

    MVector operator+(const MVector & v) const
    {
        return opVecVec(*this, v, OP_PLUS);
    }

    MVector operator-(const MVector & v) const
    {
        return opVecVec(*this, v, OP_MINUS);
    }

    MVector operator*(const MVector & v) const
    {
        return opVecVec(*this, v, OP_MULT);
    }

    MVector operator/(const MVector & v) const
    {
        return opVecVec(*this, v, OP_DIVIDE);
    }


    template <class scale_t>
    MVector operator+(const scale_t & x) const
    {
        return opVecScale(*this, x, OP_PLUS);
    }

    template <class scale_t>
    MVector operator-(const scale_t & x) const
    {
        return opVecScale(*this, x, OP_MINUS);
    }

    template <class scale_t>
    MVector operator*(const scale_t & x) const
    {
        return opVecScale(*this, x, OP_MULT);
    }

    template <class scale_t>
    MVector operator/(const scale_t & x) const
    {
        return opVecScale(*this, x, OP_DIVIDE);
    }

    template <class scale_t>
    friend MVector operator*(const scale_t & x, const MVector & v)
    {
        return v*x;
    }

    //unit operator
    MVector roundInt(void) const
    {
        MVector ret;
        for (size_t i = 0; i != m_vec.size(); ++i)
        {
            int r = (m_vec[i] > 0)? (m_vec[i] + 0.5): (m_vec[i]-0.5);
            ret.push_back(r);
        }
        return ret;
    }

    //---------control============

    template <typename SrcT>
    void attach(SrcT itr,  SrcT end)
    {
        m_vec.clear();
        while (itr != end)
            m_vec.push_back(*itr++);
    }

    //==============staic===
    static T pow(T x, int n)
    {
        assert(n >= 0);
        T ret = 1;

        for (int i = 0; i < n; ++i)
            ret *= x;
        return ret;
    }


    //================= test
    void show(void) const
    {
        using namespace std;
        typename vector_t::const_iterator itr = m_vec.begin();
        for (; itr != m_vec.end(); ++itr)
            cout << *itr << " ";

        cout << endl;
    }


    friend std::ostream& operator<<(std::ostream &os, const MVector & v)
    {
        using namespace std;
        typename vector_t::const_iterator itr = v.m_vec.begin();
        for (; itr != v.m_vec.end(); ++itr)
            os << *itr << " ";
        return os;
    }

};



#endif /* end of include guard: MVECTOR_H */
