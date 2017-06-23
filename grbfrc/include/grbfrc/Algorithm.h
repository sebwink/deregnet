#ifndef GRBFRC_ALGORITHM_H
#define GRBFRC_ALGORITHM_H

#include <vector>
#include <memory>

#include <gurobi_c++.h>

#include <grbfrc/common.h>

namespace grbfrc {

class _Algorithm {

    protected:

        GRBModel* base_model;
        FMILPObj* objective;
        std::vector<GRBVar*>* vars;
        std::vector<double>* start_solution;

        GRBCallback* callback;

    public:

        void init(GRBModel* xbase_model,
                  FMILPObj* xobjective,
                  std::vector<GRBVar*>* xvars,
                  std::vector<double>* xstart_solution);
/*
        template <typename T>
        void run(GrbfrcCallback<T>* cb);
*/
        virtual void writeSolution(FMILPSol** solution) = 0;

};

/*
template <typename T>
void _Algorithm::run(GrbfrcCallback<T>* cb) {
    // throw NotImplementedError
}
*/

}

#endif
