#include <gurobi_c++.h>
#include <grbfrc.h>

using FMILP= grbfrc::FMILP;

int main()
 {
  GRBEnv env;
  FMILP model(env);
  GRBVar x = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "x");
  GRBVar y = model.addVar(0.0, 6.0, GRB_CONTINUOUS, "y");
  model.update();
  model.addConstr(-x + y <= 4.0);
  model.addConstr(2*x + y <= 14.0);
  GRBLinExpr objNum { -2*x + y + 2};
  model.setObjNumerator(objNum);
  GRBLinExpr objDenom { x + 3*y + 4 };
  model.setObjDenominator(objDenom);
  model.setObjSense(GRB_MINIMIZE);
  model.optimize(grbfrc::Algorithm::CCT);
  model.printSolution();
 }
