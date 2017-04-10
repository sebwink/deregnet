#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <grbfrc/CharnesCooper.h>
#include <grbfrc/Dinkelbach.h>
#include <grbfrc/YGGY.h>
#include <grbfrc/FMILP.h>

namespace grbfrc
{

FMILP::FMILP(GRBEnv& env,
             GRBLinExpr objNumerator,
             GRBLinExpr objDenominator,
             int objSense)

             : baseModel { new GRBModel(env) }
 {
  objective.numerator = objNumerator;
  objective.denominator = objDenominator;
  baseModel->setObjective(objDenominator);
  objective.sense = objSense;
  objective.uniSign = UniSignStatus::maybe;
  solution = nullptr;
 }

FMILP::FMILP(GRBEnv& grbEnv)
       : FMILP(grbEnv, 0.0, 0.0, GRB_MAXIMIZE)
 { objective.uniSign = UniSignStatus::yes; }

FMILP::FMILP(GRBEnv& grbEnv,
                 GRBLinExpr objNumerator,
                 GRBLinExpr objDenominator)
       : FMILP(grbEnv, objNumerator, objDenominator, GRB_MAXIMIZE)
 { }

FMILP::FMILP(GRBEnv& grbEnv,
                 int objSense)
       : FMILP(grbEnv, 0.0, 0.0, objSense)
 { objective.uniSign = UniSignStatus::yes; }


GRBModel FMILP::baseModelCopy()
 { return *baseModel; }

FMILP FMILP::copy()
 { return *this; }

int FMILP::get(GRB_IntAttr attr)
 { return baseModel->get(attr); }

double FMILP::get(GRB_DoubleAttr attr)
 { return baseModel->get(attr); }

std::string FMILP::get(GRB_StringAttr attr)
 { return baseModel->get(attr); }

void FMILP::set(GRB_IntAttr attr, int val)
 { baseModel->set(attr, val); }

void FMILP::set(GRB_DoubleAttr attr, double val)
 { baseModel->set(attr, val); }

void FMILP::set(GRB_StringAttr attr, std::string val)
 { baseModel->set(attr, val); }

void FMILP::set(GRB_IntParam param, int val)
 { baseModel->set(param, val); }

void FMILP::set(GRB_DoubleParam param, double val)
 { baseModel->set(param, val); }

void FMILP::set(GRB_StringParam param, std::string val)
 { baseModel->set(param, val); }

void FMILP::setEnv(GRB_IntParam para, int val)
 { (baseModel->getEnv()).set(para, val); }

GRBVar& FMILP::addVar(double lowerBound,
                        double upperBound,
                        char type,
                        std::string name)
 {
  vars.push_back(new GRBVar(baseModel->addVar(lowerBound, upperBound, 0.0, type, name)));
  return *(vars.back());
 }

GRBVar& FMILP::addVar(double lowerBound,
                        double upperBound,
                        double numeratorCoeff,
                        double denominatorCoeff,
                        char type,
                        std::string name)
 {
  vars.push_back(new GRBVar(baseModel->addVar(lowerBound, upperBound, denominatorCoeff, type, name)));
  objective.numerator += numeratorCoeff * (*(vars.back()));
  objective.denominator += denominatorCoeff * (*(vars.back()));
  return *(vars.back());
 }

void FMILP::setStartSolution(GRBVar& var, double val) {
    if (startSol)
        startSol = new std::vector<double>(vars.size());
    else
        startSol->resize(vars.size());
    int index { 0 };
    for (GRBVar* varp : vars) {
        if (varp->sameAs(var))
            (*startSol)[index] = val;
        index++;
    }

}

GRBVar* FMILP::getVars()
 { return baseModel->getVars(); }

void FMILP::update()
 { baseModel->update(); }

void FMILP::setObjNumerator(GRBLinExpr& objNumerator)
 { objective.numerator = objNumerator; }

void FMILP::setObjDenominator(GRBLinExpr& objDenominator)
 { objective.denominator = objDenominator; }

int FMILP::getObjSense()
 { return objective.sense; }

void FMILP::setObjSense(int objSense)
 { objective.sense = objSense; }

void FMILP::addConstr(const GRBLinExpr lhsExpr,
                        char sense,
                        const GRBLinExpr rhsExpr,
                        std::string name)
 { baseModel->addConstr(lhsExpr, sense, rhsExpr, name); }

void FMILP::addConstr(const GRBLinExpr lhsExpr,
                        char sense,
                        GRBVar rhsVar,
                        std::string name)
 { baseModel->addConstr(lhsExpr, sense, rhsVar, name); }

void FMILP::addConstr(const GRBLinExpr lhsExpr,
                        char sense,
                        double rhsVal,
                        std::string name)
 { baseModel->addConstr(lhsExpr, sense, rhsVal, name); }

void FMILP::addConstr(GRBVar lhsVar,
                        char sense,
                        GRBVar rhsVar,
                        std::string name)
 { baseModel->addConstr(lhsVar, sense, rhsVar, name); }

void FMILP::addConstr(GRBVar lhsVar,
                        char sense,
                        double rhsVal,
                        std::string name)
 { baseModel->addConstr(lhsVar, sense, rhsVal, name); }

void FMILP::addConstr(GRBTempConstr tc,
                        std::string name)
 { baseModel->addConstr(tc, name); }

void FMILP::addConstr(const GRBLinExpr expr,
                        double lower,
                        double upper)
 {
  baseModel->addConstr(expr, GRB_GREATER_EQUAL, lower);
  baseModel->addConstr(expr, GRB_LESS_EQUAL, upper);
 }

GRBConstr* FMILP::getConstrs()
 { return baseModel->getConstrs(); }

GRBLinExpr FMILP::getRow(GRBConstr& constr)
 { return baseModel->getRow(constr); }

bool FMILP::isFLP()
 {
  if (baseModel->get(GRB_IntAttr_NumIntVars) == 0) return true;
  else return false;
 }

bool FMILP::checkUnisignance(std::string sign)
 {
  /**
  *
  * test
  */
  if (sign == "positive")
     baseModel->setObjective(objective.denominator, GRB_MINIMIZE);
  else if (sign == "negative")
     baseModel->setObjective(objective.denominator, GRB_MAXIMIZE);

  baseModel->optimize();

  if (sign == "positive" && baseModel->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
   {
    if (baseModel->get(GRB_DoubleAttr_ObjVal) > 0)
     {
      objective.uniSign = UniSignStatus::yes;
      return true;
     }
    else return false ;
   }
  else return false;
  if (sign == "negative" && baseModel->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
   {
    if (baseModel->get(GRB_DoubleAttr_ObjVal) < 0)
     {
      objective.uniSign = UniSignStatus::yes;
      return true;
     }
    else return false;
   }
  else return false;
 }

bool FMILP::isUnisignant()
 {
  /** test
  * \callgraph
  *
  *
  */
  if (checkUnisignance("positive"))
   {
    objective.uniSign = UniSignStatus::yes;
    return true;
   }
  else if (checkUnisignance("negative"))
   {
    objective.uniSign = UniSignStatus::yes;
    return true;
   }
  else
   {
    objective.uniSign = UniSignStatus::no;
    return false;
   }
 }

void FMILP::optimize(Algorithm algorithm) {
    if (algorithm == Algorithm::GCC)
        this->runGeneralizedCharnesCooper();
    else if (algorithm == Algorithm::CCT)
        this->runCharnesCooper();
    else if (algorithm == Algorithm::DTA)
        this->runDinkelbach();
    //else if (algorithm == Algorithm::OVT)
        //this->runObjVarTransform();
}

double FMILP::getObjVal() {
    return solution->objVal;
}

double FMILP::getVal(GRBVar& var) {
    int index { 0 };
    for (GRBVar* varp : vars) {
        if (varp->sameAs(var))
            return (solution->varVals)[index];
        index++;
    }
    return 0.0; // bad
}

void FMILP::runCharnesCooper(int time_limit)
 {
  baseModel->update();
  CharnesCooper chaco(*this);
  chaco.run(time_limit);
  chaco.writeSolution();
 }

void FMILP::runCharnesCooper(GRBCallback& callback,
                               int time_limit)
 {
  baseModel->update();
  CharnesCooper chaco(*this);
  chaco.run(callback, time_limit);
  chaco.writeSolution();
 }

GRBModel FMILP::getCharnesCooperTransform()
 {
  CharnesCooper chaco(*this);
  return chaco.getTransform();
 }

void FMILP::runDinkelbach(double qi,
                            int max_iter,
                            int time_limit,
                            double tol,
                            bool verbose,
                            bool logFile,
                            std::string logFileName)
 {
  baseModel->update();
  Dinkelbach dinkel(*this, qi, max_iter, time_limit, tol);
  dinkel.run(verbose, logFile, logFileName);
  dinkel.writeSolution();
 }

void FMILP::runDinkelbach(GRBCallback& callback,
                            double qi,
                            int max_iter,
                            int time_limit,
                            double tol,
                            bool verbose,
                            bool logFile,
                            std::string logFileName)
 {
  baseModel->update();
  Dinkelbach dinkel(*this, qi, max_iter, time_limit, tol);
  dinkel.run(callback, verbose, logFile, logFileName);
  dinkel.writeSolution();
 }

void FMILP::runGeneralizedCharnesCooper(int time_limit) {
    baseModel->update();
    YGGY yggy(this);
    yggy.run(time_limit);
    yggy.writeSolution();
}

void FMILP::runYGGY(GRBCallback callback,
             int time_limit)
 { }

GRBModel FMILP::getYGGYTransform()
 {
  GRBModel model(baseModel->getEnv());
  return model;
 }

void FMILP::printSolution()
 {
  if (solution)
   {
    std::cout << "\n########################### Solution ############################\n\n";
    for (int i = 0; i < baseModel->get(GRB_IntAttr_NumVars); i++)
     {
      double varval { solution->varVals[i] };
      if (varval == 0.0) varval = std::abs(varval);
      std::cout << vars[i]->get(GRB_StringAttr_VarName)
                << ": " << varval << std::endl;
     }
    std::cout << std::endl;
    std::cout << "ObjVal: " << solution->objVal << std::endl;
    std::cout << "\n#################################################################\n\n";
   }
  else
    std::cout << "There is no solution to print!" << std::endl;
 }


}       // namespace grbfrc
