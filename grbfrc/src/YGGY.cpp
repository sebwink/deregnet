#include <gurobi_c++.h>
#include <grbfrc/YGGY.h>
#include <grbfrc/Gloverizer.h>

namespace grbfrc
{

YGGY::YGGY(FMILP* fmipPtr)
             : fmip { fmipPtr },
               invalid { false },
               transformed { false },
               transformation { GRBModel((fmipPtr->baseModel)->getEnv()) } {
    fmip->update();
}

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

int YGGY::transform() {
    // denominator variable u = 1 / (dx + ey + f)
    double umax = 0.0;
    if (!getUmax(umax)) return 127;
    u = transformation.addVar(0.0, umax, 0.0, GRB_CONTINUOUS);
    // tx ~ (u*xc, xd)
    define_tx();
    transformation.update();
    FMILPObj& objective = fmip->objective;
    define_objective(objective.numerator, objective.sense);
    define_constraints(objective.denominator);
    transformation.update();
    Gloverizer gloverizer(&transformation, &u, umax);
    gloverizer.gloverize();
    return 0;
}

bool YGGY::getUmax(double& umax) {
    std::cout << "Calculating Umax ...\n " << std::endl;    // if verbose ...
    fmip->baseModel->setObjective((fmip->objective).denominator, GRB_MINIMIZE);
    fmip->baseModel->optimize();
    std::cout << "\nDone calculating Umax.\n" << std::endl;       // if verbose ...
    if (fmip->get(GRB_IntAttr_Status) != GRB_OPTIMAL) return false;
    umax = 1 / fmip->get(GRB_DoubleAttr_ObjVal);
    return true;
}

void YGGY::define_tx() {
    int vindex { 0 };
    for (auto var: fmip->vars) {
        double lb { var->get(GRB_DoubleAttr_LB) };
        double ub { var->get(GRB_DoubleAttr_UB) };
        if (var->get(GRB_CharAttr_VType) == 'C') {
            tx.push_back( new GRBVar( transformation.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS) ) );
            transformation.update();
            if (lb != -GRB_INFINITY)
                transformation.addConstr(*(tx.back()) - lb*u >= 0.0);  // get umin & use ?!
            if (ub != GRB_INFINITY)
                transformation.addConstr(*(tx.back()) - ub*u <= 0.0);  // use umax ?!
        }
        else {
            tx.push_back( new GRBVar( transformation.addVar(lb, ub, 0.0, GRB_BINARY) ) );
        }
        vindex++;
    }
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


    GRBConstr* constrs { fmip->getConstrs() };
    int numConstrs { fmip->get(GRB_IntAttr_NumConstrs) };
    for (int i = 0; i < numConstrs; i++) {
      GRBLinExpr origLhs { fmip->getRow(*constrs) };
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
        transform();
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

void YGGY::solveTransform(GRBCallback& callback) {
    if (invalid)
        printInvalidity();
    else {
        try {
              transformation.setCallback(&callback);
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

void YGGY::run(int time_limit) {
    if (invalid)
        printInvalidity();
    else {
        std::cout << "\n=========== solving FMIP via YGGY transform ===========\n\n";
        if (!transformed) transform();
        solveTransform();
        backTransformSolution();
    }
}

void YGGY::run(GRBCallback& callback, int time_limit) {
    if (invalid)
        printInvalidity();
    else {
        std::cout << "\n=========== solving FLP via Charnes-Cooper transform ===========\n\n";
        if (!transformed) transform();
        solveTransform(callback);
        backTransformSolution();
    }
}

void YGGY::writeSolution() {
    if (fmip->solution) *(fmip->solution) = solution;
    else if (!fmip->solution) fmip->solution = new FMILPSol(solution);
    else std::cout << "No solution avaiable!" << std::endl;
}

FMILPSol YGGY::getSolution() {
    return solution;
}

void YGGY::backTransformSolution() {
    int status { transformation.get(GRB_IntAttr_Status) };
    if (status == GRB_OPTIMAL) {
      solution.objVal = transformation.get(GRB_DoubleAttr_ObjVal);
      double uval { u.get(GRB_DoubleAttr_X) };
      for (int i = 0; i < fmip->get(GRB_IntAttr_NumVars); i++) {
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
    for (int i = 0; i < fmip->get(GRB_IntAttr_NumVars); i++)
        if (var.sameAs(*(fmip->vars[i]))) return i;
    return -1;
}

} // namespace grbfrc
