#ifndef GRBFRC_SOLVER_H
#define GRBFRC_SOLVER_H

#include <vector>

#include <gurobi_c++.h>

#include <grbfrc/common.h>
#include <grbfrc/Algorithm.h>
#include <grbfrc/CharnesCooper.h>
#include <grbfrc/YGGY.h>
#include <grbfrc/Dinkelbach.h>
#include <grbfrc/ObjVarTransform.h>

namespace grbfrc {

class Solver {

    private:

        GRBModel* base_model { nullptr };
        FMILPObj* objective { nullptr };
        std::vector<GRBVar*>* vars { nullptr };
        std::vector<double>* start_solution { nullptr };
        FMILPSol** solution { nullptr };

    public:

        Solver(GRBModel* base_model_ptr,
               FMILPObj* objptr,
               std::vector<GRBVar*>* varsptr,
               std::vector<double>* startptr,
               FMILPSol** solptr);

        template <typename T>
        void solve(Algorithm algorithm = Algorithm::GCC, 
                   GrbfrcCallback<T>* cb = nullptr,
                   double* objub = nullptr,
                   double* objlb = nullptr);

        template <typename T>
        void _solve(CharnesCooper* algorithm,
                    GrbfrcCallback<T>* cb = nullptr);

        template <typename T>
        void _solve(Dinkelbach* algorithm,
                    GrbfrcCallback<T>* cb = nullptr);

        template <typename T>
        void _solve(YGGY* algorithm,
                    GrbfrcCallback<T>* cb = nullptr);

        template <typename T>
        void _solve(ObjVarTransform* algorithm,
                    GrbfrcCallback<T>* cb = nullptr);

};

template <typename T>
void Solver::solve(Algorithm algorithm, GrbfrcCallback<T>* cb,
                   double* objub, double* objlb) {
    GRBEnv env = base_model->getEnv();
    if (algorithm == Algorithm::CCT) {
        CharnesCooper _algorithm = CharnesCooper(&env);
        this->_solve(&_algorithm, cb);
    }
    else if (algorithm == Algorithm::DTA) {
        Dinkelbach _algorithm = Dinkelbach();
        this->_solve(&_algorithm, cb);
    }
    else if (algorithm == Algorithm::OVT) {
        ObjVarTransform _algorithm = ObjVarTransform(&env, objub, objlb);
        this->_solve(&_algorithm, cb);
    }
    else {
        YGGY _algorithm = YGGY(&env);
        this->_solve(&_algorithm, cb);
    }

}

template <typename T>
void Solver::_solve(CharnesCooper* algorithm, GrbfrcCallback<T>* cb) {
    algorithm->init(base_model,
                    objective,
                    vars,
                    start_solution);
    algorithm->run(cb);
    // check for sanity ...
    // then:
    algorithm->writeSolution(solution);
}

template <typename T>
void Solver::_solve(Dinkelbach* algorithm, GrbfrcCallback<T>* cb) {
    algorithm->init(base_model,
                    objective,
                    vars,
                    start_solution);
    algorithm->run(cb);
    // check for sanity ...
    // then:
    algorithm->writeSolution(solution);
}

template <typename T>
void Solver::_solve(YGGY* algorithm, GrbfrcCallback<T>* cb) {
    algorithm->init(base_model,
                    objective,
                    vars,
                    start_solution);
    algorithm->run(cb);
    // check for sanity ...
    // then:
    algorithm->writeSolution(solution);
}

template <typename T>
void Solver::_solve(ObjVarTransform* algorithm, GrbfrcCallback<T>* cb) {
    algorithm->init(base_model,
                    objective,
                    vars,
                    start_solution);
    algorithm->run(cb);
    // check for sanity ...
    // then:
    algorithm->writeSolution(solution);
}

// specialized _solve methods if algorithm parameters neccessary.

}    //    namespace grbfrc

#endif    //    GRBFRC_SOLVER_H
