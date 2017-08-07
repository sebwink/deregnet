#include <gurobi_c++.h>
#include <grbfrc.h>

using FMILP = grbfrc::FMILP;

int main(void)
 {
  GRBEnv env;
  FMILP model(env);
  GRBVar x = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "x");
  GRBVar y = model.addVar(0.0, 6.0, GRB_CONTINUOUS, "y");
  model.update();
  model.addConstr(-x + y <= 4.0);
  model.addConstr(2*x + y <= 14.0);

  GRBVar u = model.addVar(0.0, 1.0, GRB_BINARY, "u");
  GRBVar v = model.addVar(0.0, 1.0, GRB_BINARY, "v");
  model.update();
  model.addConstr( u + v <= 1 );

  GRBLinExpr objNum { -2*x + y + u + 2 };
  model.setObjNumerator(objNum);
  GRBLinExpr objDenom { x + 3*y + u + v + 4 };
  model.setObjDenominator(objDenom);
  model.setObjSense(GRB_MINIMIZE);

  model.optimize(grbfrc::Algorithm::DTA);
  model.printSolution();  

  return 0;
 }
