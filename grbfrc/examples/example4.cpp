#include <gurobi_c++.h>
#include <grbfrc.h>

using FMILP = grbfrc::FMILP;

int main()
 {
  GRBEnv env;
  FMILP model(env);
  GRBVar x = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "x");
  GRBVar y = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "y");
  GRBVar z = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "z");
  GRBVar u = model.addVar(0.0, 1.0, GRB_BINARY, "u");
  GRBVar v = model.addVar(0.0, 1.0, GRB_BINARY, "v");
  model.update();
  model.addConstr( x + 3*y + z <= 10 );
  model.addConstr( 3*x + 2*y <= 7 );
  model.addConstr( u + v >= 1 );
  GRBLinExpr objNum { 2*x + y + 3*z + u + v + 6 };
  model.setObjNumerator(objNum);
  GRBLinExpr objDenom { x + 3*z + v + 4 };
  model.setObjDenominator(objDenom);
  model.setObjSense(GRB_MAXIMIZE);
  model.optimize(grbfrc::Algorithm::DTA);
  model.printSolution();
 }
