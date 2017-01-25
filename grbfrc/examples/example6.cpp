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
  model.update();
  model.addConstr( x + 3*y + z <= 10 );
  model.addConstr( 3*x + 2*y <= 7 );
  GRBLinExpr objNum { 2*x + y + 3*z + 6 };
  model.setObjNumerator(objNum);
  GRBLinExpr objDenom { x + 3*z + 4 };
  model.setObjDenominator(objDenom);
  model.setObjSense(GRB_MAXIMIZE);
  std::cout << "\nSolve with Charnes-Cooper transformation:" << std::endl;
  model.runCharnesCooper();
  model.printSolution();
  std::cout << "\nSolve with generalized Charnes-Cooper transformation:" << std::endl;
  model.runYGGY();
  model.printSolution();
  std::cout << "\nSolve with Dinkelbach algorithm: " << std::endl;
  model.runDinkelbach();
  model.printSolution();
 }
