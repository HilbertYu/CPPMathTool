#ifndef MVECTORNESTLOOP_H
#define MVECTORNESTLOOP_H

#include <vector>
#include "MVector.h"

template <typename scale_t>
class MVectorNestLoop
{
public:
    class LoopCtx
    {
    public:
        LoopCtx(scale_t b = 0, scale_t e = 0, scale_t sp = 1):
            begin(b),
            end(e),
            step(sp)
        {

        }

        scale_t begin;
        scale_t end;
        scale_t step;
    };

    MVectorNestLoop(): m_bds() { }
    MVectorNestLoop(const std::vector<LoopCtx> &bds): m_bds(bds) { }

    void addLoopLevel(scale_t b, scale_t e, scale_t sp = 1)
    {
        m_bds.push_back(LoopCtx(b,e,sp));
    }

    void clean(void)
    {
        m_bds.clear();
    }

    template <typename Func_t>
    void doLoop(Func_t func) const
    {
        MVector<scale_t> x;
        x.zeros(m_bds.size());

        ImpDoLoop(m_bds.size() - 1, x, func);
    }

private:
    std::vector<LoopCtx> m_bds;

    template <typename Func_t>
    void ImpDoLoop(int loop_lv, MVector<scale_t> & itr_x, Func_t & func) const
    {
        if (loop_lv < 0)
            return;

        const std::vector<LoopCtx> & bds = m_bds;

        for (scale_t v = bds[loop_lv].begin; v < bds[loop_lv].end; v = v + bds[loop_lv].step)
        {
            itr_x[loop_lv] = v;

            if (loop_lv > 0)
            {
                ImpDoLoop(loop_lv - 1, itr_x, func);
            }

            if (loop_lv == 0)
            {
                func(itr_x);
            }
        }

    }

};

#endif /* end of include guard: MVECTORNESTLOOP_H */
