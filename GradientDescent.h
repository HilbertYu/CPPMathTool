#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H
#include "MFunction.h"

template <size_t dim, typename scale_t>
class GradientDescent
{
    const int m_max_itr_times;
    const double m_zero_tol;
    const double m_end_tol;

    MVector<double> m_cur_ret_x;
    int m_cur_itr_times;
public:
    GradientDescent():
        m_max_itr_times(20),
        m_zero_tol(0.5),
        m_end_tol(0.5),
        m_cur_ret_x(),
        m_cur_itr_times(0)
    {

    }

    const MVector<double> & getRetX(void) const
    {
        return m_cur_ret_x;
    }

    int getRetItrTimes(void) const
    {
        return m_cur_itr_times;
    }

    virtual ~GradientDescent(void) {}

    MVector<double> gradent_descent_scale_get_next_step(
            const MFunction<dim, scale_t> & mf,
            const MVector<double> & x,
            const MVector<double> & z)
    {
        double g1 = mf.func(x);
        double extream_val = g1;
        int da = 16;

        const MVector<double> & rx = (x - da*z).roundInt();
        double tmp_val = mf.func(rx);

        while (g1 <= tmp_val && da > 0)
        {
            da /= 2;
            const MVector<double> & rx = (x - da*z).roundInt();
            tmp_val = mf.func(rx);
        }

        MVector<double> next_x  = (x - da*z).roundInt();

        if (da == 0)
        {
            return x;
        }

        for (int r = 0; r < da; ++r)
        {
            const MVector<double> & rx = (x - r*z).roundInt();
            double tmp_val = mf.func(rx);

            if (tmp_val < extream_val)
            {
                extream_val = tmp_val;
                next_x = rx;
            }
        }

        return next_x;
    }

    void run(const MFunction<dim, scale_t> & mf, const MVector<scale_t> & init_x)
    {

        MVector<double> x = init_x;
        MVector<double> x_ret;

        double g1 = 0;
        double z0 = 0;

        MVector<double> z;

        int itr_times = 0;
        for (int i = 0; i < m_max_itr_times; ++i, ++itr_times)
        {
#ifdef DEBUG_GD
            printf("--- frame %d/%d----\n", i, N);
#endif
            {
                for (size_t r = 0; r < x.size(); ++r)
                {
                    if (abs(x[r]) > 128)
                    {
                        itr_times = -1;
                        break;
                    }
                }
            }

            g1 = mf.func(x);
            z  = mf.gradient(x);
            z0 = z.norm();

#ifdef DEBUG_GD
            cout << "g1 = " << g1 << endl;
            cout << z << endl;
            cout << "grad norm = " << z0 << endl;
#endif

            if (z0 < m_zero_tol)
            {
#ifdef DEBUG_GD
                cout << "norm(gd) < tol: " << z0 << endl;
#endif
                x_ret = x;
                break;
            }

            //normalize
            z = z/z0;
#ifdef DEBUG_GD
            cout << "grad (normalize) = " << z0 << endl;
            cout << z << endl;
#endif

            MVector<double> next_x =
                gradent_descent_scale_get_next_step(mf, x, z);

#ifdef DEBUG_GD
            cout << "next x = ";
            cout << next_x << endl;
#endif

            double next_val = mf.func(next_x);
            double err = abs(next_val - g1);

#ifdef DEBUG_GD
            cout << "err = " << err << endl;
#endif
            if (err < m_end_tol)
            {
#ifdef DEBUG_GD
                cout << "END" << endl;
#endif
                x_ret = x;
                break;
            }
            x = next_x;
        }


        if (itr_times > m_max_itr_times - 2)
        {
            printf("over itr times = %d\n", itr_times);
            exit(0);
        }

        m_cur_ret_x = x_ret;
        m_cur_itr_times = itr_times;

    }

private:
    GradientDescent(const GradientDescent & rhs);
    GradientDescent & operator = (const GradientDescent & rhs);
        /* data_m */
}; /* end of class GradientDescent */


#endif /* end of include guard: GRADIENTDESCENT_H */
