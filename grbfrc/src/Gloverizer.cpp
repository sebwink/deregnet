#include <gurobi_c++.h>

#include <grbfrc/Gloverizer.h>


namespace grbfrc {

Gloverizer::Gloverizer(GRBModel* xmodel, GRBVar* xu, double umax) 
                       : model { xmodel },
					     u { xu }, 
						 U { umax } {
    // variables
    GRBVar* var { model->getVars() };
    double L { u->get(GRB_DoubleAttr_LB) };
    int numVars { model->get(GRB_IntAttr_NumVars) };
    for (int i = 0; i < numVars; i++) {
        if (var->get(GRB_CharAttr_VType) != 'C') {
            glvrvr[var] = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, 'C');
            binvars.push_back(var);
            model->update();
            model->addConstr( glvrvr[var] <= U*(*var) );
            model->addConstr( glvrvr[var] >= L*(*var) );
            model->addConstr( *u - U*(1 - (*var)) <= glvrvr[var] );
            model->addConstr( *u - L*(1 - (*var)) >= glvrvr[var] );
        }
        var++;
    }
    // constraints
    GRBQConstr* qconstr { model->getQConstrs() };
    for (int q = 0; q < model->get(GRB_IntAttr_NumQConstrs); q++) {
        qconstrs.push_back(qconstr);
        GRBQuadExpr qexpr { model->getQCRow(*qconstr) };
        qconstr2newlhs[qconstr] = qexpr.getLinExpr();
        qconstr++;
    }
    // objects
    GRBQuadExpr oldobj { model->getObjective() };
    newobj = oldobj.getLinExpr();
}

void Gloverizer::gloverize() {
    for (auto var : binvars) {
        for (auto qconstr : qconstrs)
            handle_qconstr(var, qconstr);
        handle_objective(var);
    }
    for (auto qconstr : qconstrs) {
        model->addConstr(qconstr2newlhs[qconstr],
                         qconstr->get(GRB_CharAttr_QCSense),
                         qconstr->get(GRB_DoubleAttr_QCRHS));
        model->remove(*qconstr);
    }
    model->setObjective(newobj, model->get(GRB_IntAttr_ModelSense));
    model->update();
}

void Gloverizer::handle_qconstr(GRBVar* var, GRBQConstr* qconstr) {
    GRBQuadExpr qexpr { model->getQCRow(*qconstr) };
    for (unsigned int i = 0; i < qexpr.size(); i++) {
        if (qexpr.getVar1(i).sameAs(*var) || qexpr.getVar2(i).sameAs(*var))
            qconstr2newlhs[qconstr] += qexpr.getCoeff(i) * glvrvr[var];
    }
}

void Gloverizer::handle_objective(GRBVar* var) {
    GRBQuadExpr oldobj { model->getObjective() };
    for (unsigned int i = 0; i < oldobj.size(); i++) {
        if (oldobj.getVar1(i).sameAs(*var) || oldobj.getVar2(i).sameAs(*var))
            newobj += oldobj.getCoeff(i) * glvrvr[var];
    }
}


} // namespace grbfrc
