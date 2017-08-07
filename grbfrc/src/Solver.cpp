
#include <grbfrc/Solver.h>

namespace grbfrc {

Solver::Solver(GRBModel* base_model_ptr,
               FMILPObj* objptr,
               std::vector<GRBVar*>* varsptr,
               std::vector<double>* startptr,
               FMILPSol** solptr)

    : base_model { base_model_ptr },
      objective { objptr },
      vars { varsptr },
      start_solution { startptr },
      solution { solptr }

   { }

}
