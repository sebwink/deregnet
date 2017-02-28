// --------------------------------------------------------------------------
//                grbfrc -- Mixed-integer fractional programming
// --------------------------------------------------------------------------
// Copyright Sebastian Winkler --- Eberhard Karls University Tuebingen, 2016
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Sebastian Winkler $
// $Authors: Sebastian Winkler $
// --------------------------------------------------------------------------
//

#ifndef FMILP_H
#define FMILP_H

#include <string>
#include <vector>

#include <gurobi_c++.h>

namespace grbfrc
{

class CharnesCooper;
class Dinkelbach;
class YGGY;


enum Algorithm {
    CCT, /*< Charnes-Cooper transformation */
    GCC, /*< Generalized Charnes-Cooper transform */
    DTA, /*< Dinkelbach-type algorithm */
    OVT  /*< Objective-Variable transform */
};

// UniSignStatus ===================================================================== //
enum UniSignStatus
 {
  maybe = 0,
  yes = 1,
  no = -1
 };
// FMILPObj ======================================================================== //
struct FMILPObj
 {
  GRBLinExpr numerator;
  GRBLinExpr denominator;
  int sense;
  UniSignStatus uniSign;
 };
// FMILPSol ======================================================================== //
struct FMILPSol
 {
  double objVal;
  std::vector<double> varVals;
 };
// FMILP =========================================================================== //
class FMILP
 {

    friend class CharnesCooper;
    friend class Dinkelbach;
    friend class YGGY;
    friend class Gloverizer;

  private:

    GRBModel* baseModel;
    FMILPObj objective;
    FMILPSol* solution = nullptr;
    std::vector<GRBVar*> vars;

  public:

    // =============================================================================== //
    // FMILP base ctor                                                               //
    // =============================================================================== //
    FMILP(GRBEnv& env,
            GRBLinExpr objNumerator,
            GRBLinExpr objDenominator,
            int objSense);
    // derived FMILP ctors ========================================================= //
    FMILP(GRBEnv& grbEnv);
    // =============================================================================== //
    FMILP(GRBEnv& grbEnv,
            GRBLinExpr objNumerator,
            GRBLinExpr objDenominator);
    // =============================================================================== //
    FMILP(GRBEnv& grbEnv,
            int objSense);
    // =============================================================================== //
    // FMILP dtor                                                                    //
    // =============================================================================== //
    // ~FMILP();
    // =============================================================================== //
    //       methods                                                                   //
    // =============================================================================== //
    // get copies of base model or the FMILP itself ================================ //
    GRBModel baseModelCopy();
    FMILP copy();
    // get attributes of the base model ============================================== //
    int get(GRB_IntAttr attr);
    double get(GRB_DoubleAttr attr);
    std::string get(GRB_StringAttr attr);
    // set attributes of the base model ============================================== //
    void set(GRB_IntAttr attr, int val);
    void set(GRB_DoubleAttr attr, double val);
    void set(GRB_StringAttr attr, std::string val);
    // set parameters of the base models ============================================= //
    void set(GRB_IntParam attr, int val);
    void set(GRB_DoubleParam attr, double val);
    void set(GRB_StringParam attr, std::string val);
    // set parameters of the base model's environment ================================ //
    void setEnv(GRB_IntParam para, int val);
    // add variables to the model ==================================================== //
    GRBVar& addVar(double lowerBound,
                   double upperBound,
                   char type,
                   std::string name = "");
    GRBVar& addVar(double lowerBound,
                   double upperBound,
                   double numeratorCoeff,
                   double denominatorCoeff,
                   char type,
                   std::string name = "");
    // get variables of the model ==================================================== //
    GRBVar* getVars();
    // update model (necessary before defining constraints with newly added variables) //
    void update();
    // Defining the objective function =============================================== //
    void setObjNumerator(GRBLinExpr& objNumerator);
    void setObjDenominator(GRBLinExpr& objDenominator);
    int getObjSense();
    void setObjSense(int obj_sense);
    // add constraints to the model ================================================== //
    void addConstr(const GRBLinExpr lhsExpr,
                   char sense,
                   const GRBLinExpr rhsExpr,
                   std::string name = "");
    void addConstr(const GRBLinExpr lhsExpr,
                   char sense,
                   GRBVar rhsVar,
                   std::string name = "");
    void addConstr(const GRBLinExpr lhsExpr,
                   char sense,
                   double rhsVal,
                   std::string name = "");
    void addConstr(GRBVar lhsVar,
                   char sense,
                   GRBVar rhsVar,
                   std::string name = "");
    void addConstr(GRBVar lhsVar,
                   char sense,
                   double rhsVal,
                   std::string name = "");
    void addConstr(GRBTempConstr tc,
                   std::string name = "");
    void addConstr(const GRBLinExpr expr,
                   double lower,
                   double upper);
    // get constraints of the model ================================================== //
    GRBConstr* getConstrs();
    GRBLinExpr getRow(GRBConstr& constr);
    // check if model is a pure continuous problem =================================== //
    bool isFLP();
    // check if denominator has unfirm sign in feasible region ======================= //
    bool checkUnisignance(std::string sign);
    bool isUnisignant();

    //
    void optimize(Algorithm algorithm = Algorithm::GCC);
    double getObjVal();
    double getVal(GRBVar& var);
    // =============================================================================== //
    // algorithms
    // =============================================================================== //
    // Charnes-Cooper transformation                                                   //
    // =============================================================================== //
    // A. Charnes, W.W. Cooper: Programming with linear fractional functionals         //
    //                          Nav. Res. Log. Q, 1962 (9) 3-4 (181-186)               //
    // =============================================================================== //
    // run Charnes-Cooper algorithm // =============================================== //
    void runCharnesCooper(int time_limit = 1200);
    void runCharnesCooper(GRBCallback& callback,
                          int time_limit = 1200);
    // get Charnes-Cooper transform // =============================================== //
    GRBModel getCharnesCooperTransform();
    // =============================================================================== //
    // Dinkelbach's algorithm                                                          //
    // =============================================================================== //
    // run Dinkelbach's algorithm // ================================================= //
    void runDinkelbach(double qi = 0.0,
                       int max_iter = 10,
                       int time_limit = 1200,
                       double tol = 0.01,
                       bool verbose = true,
                       bool logFile = true,
                       std::string logFileName = "grbfrc.log");
    void runDinkelbach(GRBCallback& callback,
                       double qi = 0.0,
                       int max_iter = 10,
                       int time_limit = 1200,
                       double tol = 0.01,
                       bool verbose = true,
                       bool logFile = true,
                       std::string logFileName = "grbfrc.log"); 
    // =============================================================================== //
    // YGGY transformation                                                             //
    // =============================================================================== //
    // D. Yue, G. Guillen-Gosalbez, F. You:                                            //
    // Global Optimization of Large-Scale Mixed Integer Linear Fractional              //
    // Programming Problems: A Reformulation-Linearization Method and Process          //
    // Scheduling Applications, AIChE J., 2013 (59) 11 (4255-4272)                     //
    // =============================================================================== //
    // run YGGY algorithm //
    void runGeneralizedCharnesCooper(int time_limit = 1200);
    void runYGGY(GRBCallback callback,
                 int time_limit = 1200);
    GRBModel getYGGYTransform();
    // =============================================================================== //
    // print solution (intented for small models only!)
    void printSolution();

 };     // class FMILP

}       // namespace grbfrc

#endif  // FMILP_H
