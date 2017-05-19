#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <gurobi_c++.h>

#include <grbfrc/Dinkelbach.h>

namespace grbfrc {

Dinkelbach::Dinkelbach(double qi,
                       int maxit,
                       int maxsec,
                       double tolerance)
                     : q_init { qi },
                       max_iter { maxit },
                       time_limit { maxsec },
                       tol { tolerance }
                    {  }

void Dinkelbach::set_q_init(double qi)
 { q_init = qi; }

void Dinkelbach::set_max_iter(int max)
 { max_iter = max; }

void Dinkelbach::set_time_limit(int seconds)
 { time_limit = seconds; }

void Dinkelbach::set_tolerance(double tolerance)
 { tol = tolerance; }

void Dinkelbach::_run(bool verbose,
                      bool logFile,
                      std::string logFileName)

 {

  // set verbosity and logging behaviour ============================================== //
  if (verbose)
   { std::cout << "\n============ Solving FMILP with Dinkelbach's algorithm ============\n"; }
  else
   { (base_model->getEnv()).set(GRB_IntParam_LogToConsole, 0); }

  if (logFile)
   { (base_model->getEnv()).set(GRB_StringParam_LogFile, logFileName); }
  else
   { (base_model->getEnv()).set(GRB_StringParam_LogFile, ""); }

  // get model data =================================================================== //
  GRBLinExpr& objNumerator = objective->numerator;
  GRBLinExpr& objDenominator = objective->denominator;
  int sense = objective->sense;

  // algorithm parameters ============================================================= //
  int iter = 0;
  double q { q_init };
  double Fq;
  bool done { false };

  // solve Dinkelbach iterations ====================================================== //
  double sign { 1.0 };
  if (sense == GRB_MINIMIZE) { sign = -1.0; }
  while (!done && iter < max_iter)
   {
    iter++;
    if (verbose) { printIteration(iter); }
    GRBLinExpr currentObj = sign * objNumerator - q * objDenominator;
    base_model->setObjective(currentObj, GRB_MAXIMIZE);
    base_model->optimize();
    int status = base_model->get(GRB_IntAttr_Status);
    if (status != GRB_OPTIMAL) { printProblem(status, iter); }
    else
     {
      Fq = base_model->get(GRB_DoubleAttr_ObjVal);
      double numVal { sign * objNumerator.getValue() };
      double denomVal { objDenominator.getValue() };
      q = numVal / denomVal;
      std::cout << "\n     q = " << q
                << " ,  |F(q)| = " << std::abs(Fq) << std::endl;
      if (std::abs(Fq) < tol)
       {
        done = true;
        if (sense == GRB_MINIMIZE) { q = -q; }
        // register solution
        solution.objVal = q;
        for (unsigned int i = 0; i < vars->size(); i++)
          solution.varVals.push_back( ((*vars)[i])->get(GRB_DoubleAttr_X) );
       }
     }
   }
 }


void Dinkelbach::writeSolution(FMILPSol** xsolution)
 { 
  if (*xsolution) **xsolution = solution;
  else if (!*xsolution) *xsolution = new FMILPSol(solution);
  else std::cout << "No solution avaiable!" << std::endl;
 }

FMILPSol Dinkelbach::getSolution()
 { return solution; }

void Dinkelbach::printProblem(int status,
                              int iter)
     {
      std::cout << "\nModel iteration " << iter
                << " is not solvable, GRB_STATUS: " << status << "\n";
     }

void Dinkelbach::printIteration(int iter)
     {
      std::cout << "\n==================== Dinkelbach iteration "
                << iter << " ===================\n\n";
     }

}  // namespace grbfrc
