#include <grbfrc/Algorithm.h>

namespace grbfrc  {

void _Algorithm::init(GRBModel* xbase_model,
                      FMILPObj* xobjective,
                      std::vector<GRBVar*>* xvars,
                      std::vector<double>* xstart_solution) {
    base_model = xbase_model;
    objective = xobjective;
    vars = xvars;
    start_solution = xstart_solution;
}
        

}
