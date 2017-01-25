#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <gurobi_c++.h>
#include <grbfrc/Dinkelbach.h>

namespace grbfrc
{

Dinkelbach::Dinkelbach(FMILP& fmipRef,
                       double qi,
                       int maxit,
                       int maxsec,
                       double tolerance)
                     : fmip { &fmipRef },
                       q_init { qi },
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

void Dinkelbach::run(bool verbose,
                     bool logFile,
                     std::string logFileName)

 {

  // set verbosity and logging behaviour ============================================== //
  if (verbose)
   { std::cout << "\n============ Solving FMILP with Dinkelbach's algorithm ============\n"; }
  else
   { ((fmip->baseModel)->getEnv()).set(GRB_IntParam_LogToConsole, 0); }

  if (logFile)
   { ((fmip->baseModel)->getEnv()).set(GRB_StringParam_LogFile, logFileName); }
  else
   { ((fmip->baseModel)->getEnv()).set(GRB_StringParam_LogFile, ""); }

  // get model data =================================================================== //
  GRBModel* model = fmip->baseModel;
  FMILPObj& objective = fmip->objective;
  GRBLinExpr& objNumerator = objective.numerator;
  GRBLinExpr& objDenominator = objective.denominator;
  int sense = objective.sense;

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
    model->setObjective(currentObj, GRB_MAXIMIZE);
    model->optimize();
    int status = model->get(GRB_IntAttr_Status);
    if (status != GRB_OPTIMAL) { printProblem(status, iter); }
    else
     {
      Fq = model->get(GRB_DoubleAttr_ObjVal);
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
        for (unsigned int i = 0; i < (fmip->vars).size(); i++)
          solution.varVals.push_back( ((fmip->vars)[i])->get(GRB_DoubleAttr_X) );
       }
     }
   }
 }


void Dinkelbach::run(GRBCallback& callback,
                     bool verbose,
                     bool logFile,
                     std::string logFileName)

 {

  // set verbosity and logging behaviour ============================================== //
  if (verbose)
   { std::cout << "\n============ Solving FMILP with Dinkelbach's algorithm ============\n"; }
  else
   { ((fmip->baseModel)->getEnv()).set(GRB_IntParam_LogToConsole, 0); }

  if (logFile)
   { ((fmip->baseModel)->getEnv()).set(GRB_StringParam_LogFile, logFileName); }
  else
   { ((fmip->baseModel)->getEnv()).set(GRB_StringParam_LogFile, ""); }

  // get model data =================================================================== //
  GRBModel* model = fmip->baseModel;
  FMILPObj& objective = fmip->objective;
  GRBLinExpr& objNumerator = objective.numerator;
  GRBLinExpr& objDenominator = objective.denominator;
  int sense = objective.sense;

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
    model->setObjective(currentObj, GRB_MAXIMIZE);
    model->setCallback(&callback);
    model->optimize();
    int status = model->get(GRB_IntAttr_Status);
    if (status != GRB_OPTIMAL) { printProblem(status, iter); }
    else
     {
      Fq = model->get(GRB_DoubleAttr_ObjVal);
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
        for (unsigned int i = 0; i < (fmip->vars).size(); i++)
          solution.varVals.push_back( ((fmip->vars)[i])->get(GRB_DoubleAttr_X) );
       }
     }
   }
 }


void Dinkelbach::writeSolution()
 { 
  if (fmip->solution) *(fmip->solution) = solution;
  else if (!fmip->solution) fmip->solution = new FMILPSol(solution);
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
