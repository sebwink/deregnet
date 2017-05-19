#ifndef GRBFRC_COMMON_H
#define GRBFRC_COMMON_H

namespace grbfrc {

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

}    //    namespace grbfrc

#endif    //    GRBFRC_COMMON_H
