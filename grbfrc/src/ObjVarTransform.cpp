#include <gurobi_c++.h>

#include <grbfrc/ObjVarTransform.h>

namespace grbfrc {

ObjVarTransform::ObjVarTransform(GRBEnv* env, double* objub, double* objlb)

                                 : invalid { false },
                                   transformed { false },
                                   transformation { GRBModel(*env) },
                                   umax { objub },
                                   umin { objlb } { }

void ObjVarTransform::printInvalidity() {
   // if (no umax available) "no umax"
   // if (no umin available) "no umin"
   // else if (continuous vars in denominator) "continuous variables in denominator"
}

void ObjVarTransform::copy_model() {

    for (int i = 0; i < base_model->get(GRB_IntAttr_NumVars); ++i) {
        GRBVar* varp = vars->operator[](i);
        tx.push_back(new GRBVar( transformation.addVar(varp->get(GRB_DoubleAttr_LB),
                                                       varp->get(GRB_DoubleAttr_UB),
                                                       0.0,
                                                       varp->get(GRB_CharAttr_VType)) ));
    }

    GRBConstr* constrs { base_model->getConstrs() };
    int numConstrs { base_model->get(GRB_IntAttr_NumConstrs) };
    for (int k = 0; k < numConstrs; k++) {
      GRBLinExpr origLhs { base_model->getRow(*constrs) };
      double b { constrs->get(GRB_DoubleAttr_RHS) };
      GRBLinExpr lhs = -b;
      for (unsigned int i = 0; i < origLhs.size(); i++) {
        GRBVar var = origLhs.getVar(i);
        int j = getIndex(var);
        double coeff = origLhs.getCoeff(i);
        lhs += coeff * ( *tx[j] );
       }
      char constrSense { constrs->get(GRB_CharAttr_Sense) };
      transformation.addConstr(lhs, constrSense, 0.0);
      constrs++;
     }
}

void ObjVarTransform::add_u_constraint() {
    const GRBLinExpr& objnum { objective->numerator };
    GRBLinExpr transformed_objnum { objnum.getConstant() };
    for (unsigned int i = 0; i < objnum.size(); ++i) {
        GRBVar var = objnum.getVar(i);
        transformed_objnum += objnum.getCoeff(i) * (*tx[getIndex(var)]);
    }
    const GRBLinExpr& objdenom { objective->denominator };
    GRBLinExpr transformed_objdenom { objdenom.getConstant() };
    for (unsigned int i = 0; i < objdenom.size(); ++i) {
        GRBVar glover_var = transformation.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, 'C');
        transformation.update();
        GRBVar var = objdenom.getVar(i);
        GRBVar* varp = tx[getIndex(var)];
        double U = *umax;
        double L = *umin;
        transformation.addConstr( glover_var <= U*(*varp) );
        transformation.addConstr( glover_var >= L*(*varp) );
        transformation.addConstr( u - U*(1 - (*varp)) <= glover_var );
        transformation.addConstr( u - L*(1 - (*varp)) >= glover_var );

        transformed_objdenom += objdenom.getCoeff(i) * glover_var;

    }
    transformation.addConstr(transformed_objdenom <= transformed_objnum);
}

GRBModel ObjVarTransform::getTransform() {
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

void ObjVarTransform::solveTransform() {
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

void ObjVarTransform::backTransformSolution() {
    int status { transformation.get(GRB_IntAttr_Status) };
    if (status == GRB_OPTIMAL) {
      solution.objVal = transformation.get(GRB_DoubleAttr_ObjVal);
      for (int i = 0; i < base_model->get(GRB_IntAttr_NumVars); i++)
          solution.varVals.push_back( tx[i]->get(GRB_DoubleAttr_X) );
    }
    else std::cout << "GRB_Status : " << transformation.get(GRB_IntAttr_Status) << std::endl;
}

void ObjVarTransform::writeSolution(FMILPSol** xsolution) {
    if (*xsolution) **xsolution = solution;
    else if (!*xsolution) *xsolution = new FMILPSol(solution);
    else std::cout << "No solution avaiable!" << std::endl;
}

int ObjVarTransform::getIndex(GRBVar& var) {
    for (int i = 0; i < base_model->get(GRB_IntAttr_NumVars); i++)
        if (var.sameAs(*((*vars)[i]))) return i;
    return -1;
}

}         // namespace grbfrc
