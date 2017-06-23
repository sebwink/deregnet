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

#include <memory>

#include <gurobi_c++.h>
#include <grbfrc/YGGY.h>
#include <grbfrc/Gloverizer.h>

namespace grbfrc
{

YGGY::YGGY(GRBEnv* env) : transformation { GRBModel(*env) } { }

void YGGY::printInvalidity() {
   std::cout << "\n!!! === this->printInvalidity() doesn't make sense' === !!!\n\n";
}

/*
 *  max (cc*xc + cd*xd + d)/(ec*xc + ed*xd + f)
 *
 *  s.t.
 *
 *      Ac*xc + Ad*xd = b
 *      xc in R
 *      xd in Z
 *
 */

/*
 *  u = 1 / (ec*xc + ed*xd + f)
 *  z = u*xc = xc / (ec*xc + ed*xd + f)
 *
 */

/*
 *   max (cc*z + cd*(u*xd) + d*u)
 *
 *   s.t.
 *
 *       Ac*z + Ad*(u*xd) - b*u = 0
 *       ec*z + ed*(u*xd) + f*u = 1
 *       z, u in R
 *       xd in Z
 *
 */

// Gloverize( u * xd ) --> w = u * xd

/*
 *   max (cc*z + cd*w + d*u)
 *
 *   s.t.
 *
 *       Ac*z + Ad*w - b*u = 0
 *       ec*z + ed*w + f*u = 1
 *       GloverConstraints(w, xd, u)
 *       u,w,z in R
 *       xd in Z
 *
 */

bool YGGY::getUmax(double& umax) {
    std::cout << "Calculating Umax ...\n " << std::endl;    // if verbose ...
    base_model->setObjective(objective->denominator, GRB_MINIMIZE);
    base_model->optimize();
    std::cout << "\nDone calculating Umax.\n" << std::endl;       // if verbose ...
    if (base_model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) return false;
    umax = 1 / base_model->get(GRB_DoubleAttr_ObjVal);
    return true;
}

void YGGY::define_tx() {
    int vindex { 0 };
    for (auto var: *vars) {
        double lb { var->get(GRB_DoubleAttr_LB) };
        double ub { var->get(GRB_DoubleAttr_UB) };
        if (var->get(GRB_CharAttr_VType) == 'C') {
            tx.push_back( new GRBVar( transformation.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS) ) );
            transformation.update();
            if (lb != -GRB_INFINITY)
                transformation.addConstr(*(tx.back()) - lb*u >= 0.0);  // get umin & use ?!
            if (ub != GRB_INFINITY)
                transformation.addConstr(*(tx.back()) - ub*u <= 0.0);  // use umax ?!
            // set start solution ... z_start = u_start * var_start
        }
        else {
            tx.push_back( new GRBVar( transformation.addVar(lb, ub, 0.0, GRB_BINARY) ) );
            if (start_solution) {
                tx.back()->set(GRB_DoubleAttr_Start, (*start_solution)[vindex]);
            }
        }
        vindex++;
    }
}

void YGGY::do_the_rest(double umax) {
    std::cout << "Transform objective ..." << std::endl;
    define_objective(objective->numerator, objective->sense);
    std::cout << "Transform constraints ... " << std::endl;
    define_constraints(objective->denominator);
    transformation.update();
    std::cout << "Gloverize..." << std::endl;
    Gloverizer gloverizer(&transformation, &u, umax);
    gloverizer.gloverize();
    std::cout << "Done with transforming..." << std::endl;
}

void YGGY::define_objective(GRBLinExpr& objNumerator, int objSense) {
    double d { objNumerator.getConstant() };
    GRBQuadExpr transformed_objective { d*u };
    for (unsigned int i = 0; i < objNumerator.size(); i++) {
      GRBVar var = objNumerator.getVar(i);
      int j = getIndex(var);
      double coeff = objNumerator.getCoeff(i);
      if (var.get(GRB_CharAttr_VType) == 'C')
          transformed_objective += coeff * ( *tx[j] );
      else
          transformed_objective += coeff * u * ( *tx[j] );
     }
    transformation.setObjective(transformed_objective, objSense);
}

void YGGY::define_constraints(GRBLinExpr& objDenominator) {
    double f { objDenominator.getConstant() };
    GRBQuadExpr lhs { f * u };
    for (unsigned int i = 0; i < objDenominator.size(); i++) {
      GRBVar var = objDenominator.getVar(i);
      int j = getIndex(var);
      double coeff = objDenominator.getCoeff(i);
      if (var.get(GRB_CharAttr_VType) == 'C')
          lhs += coeff * ( *tx[j] );
      else
          lhs += coeff * u * ( *tx[j] );
    }
    transformation.addQConstr(lhs, GRB_EQUAL, 1.0);


    GRBConstr* constrs { base_model->getConstrs() };
    int numConstrs { base_model->get(GRB_IntAttr_NumConstrs) };
    for (int k = 0; k < numConstrs; k++) {
      GRBLinExpr origLhs { base_model->getRow(*constrs) };
      double b { constrs->get(GRB_DoubleAttr_RHS) };
      lhs = - b * u ;
      for (unsigned int i = 0; i < origLhs.size(); i++) {
        GRBVar var = origLhs.getVar(i);
        int j = getIndex(var);
        double coeff = origLhs.getCoeff(i);
        if (var.get(GRB_CharAttr_VType) == 'C')
            lhs += coeff * ( *tx[j] );
        else
            lhs += coeff * u * ( *tx[j] );
       }
      char constrSense { constrs->get(GRB_CharAttr_Sense) };
      transformation.addQConstr(lhs, constrSense, 0.0);
      constrs++;
     }
    transformation.update();
}

GRBModel YGGY::getTransform() {
    if (invalid) {
        printInvalidity();
        return transformation;
    }
    else if (transformed)
        return transformation;
    else {
        //transform();
        return transformation;
    }
}

void YGGY::solveTransform() {
    if (invalid)
        printInvalidity();
    else {
        try {
            transformation.optimize();
            std::cout << std::endl;
        }
        catch (GRBException e) {
            std::cout << "Gurobi error: " << e.getMessage() << "\n\n";
        }
        catch (...) {
            std::cout << "Error while attempting to optimize transform ... \n\n";
        }
    }
}

void YGGY::writeSolution(FMILPSol** xsolution) {
    if (*xsolution) **xsolution = solution;
    else if (!*xsolution) *xsolution = new FMILPSol(solution);
    else std::cout << "No solution avaiable!" << std::endl;
}

/*
FMILPSol YGGY::getSolution() {
    return solution;
}
*/

void YGGY::backTransformSolution() {
    int status { transformation.get(GRB_IntAttr_Status) };
    if (status == GRB_OPTIMAL) {
      solution.objVal = transformation.get(GRB_DoubleAttr_ObjVal);
      double uval { u.get(GRB_DoubleAttr_X) };
      for (int i = 0; i < base_model->get(GRB_IntAttr_NumVars); i++) {
        double value { tx[i]->get(GRB_DoubleAttr_X) };
        if (tx[i]->get(GRB_CharAttr_VType) == 'C')
            solution.varVals.push_back( value / uval );
        else
            solution.varVals.push_back( value );
      }
    }
    else std::cout << "GRB_Status : " << transformation.get(GRB_IntAttr_Status) << std::endl;
}

int YGGY::getIndex(GRBVar& var) {
    for (int i = 0; i < base_model->get(GRB_IntAttr_NumVars); i++)
        if (var.sameAs(*((*vars)[i]))) return i;
    return -1;
}

} // namespace grbfrc
