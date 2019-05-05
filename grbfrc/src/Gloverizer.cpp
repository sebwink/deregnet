#include <gurobi_c++.h>

#include <grbfrc/Gloverizer.h>


namespace grbfrc {

Gloverizer::Gloverizer(GRBModel* xmodel, GRBVar* xu, double umax)
                       : model { xmodel },
					     u { xu }, 
                         U { umax } { }

void Gloverizer::gloverize() {
    add_glover_variables();
    std::cout << "Linearizing quadratic constraints and objective..." << std::endl;
    GRBQConstr* qconstr { model->getQConstrs() };
    for (int q = 0; q < model->get(GRB_IntAttr_NumQConstrs); q++) {
        linearize_qconstr(qconstr);
        qconstr++;
    }
    linearize_objective();
    model->update();
}

void Gloverizer::add_glover_variables() {
    std::cout << "Adding Glover variables ..." << std::endl;
    GRBVar* var { model->getVars() };
    double L { u->get(GRB_DoubleAttr_LB) };
    int numVars { model->get(GRB_IntAttr_NumVars) };
    for (int i = 0; i < numVars; i++) {
        if (var->get(GRB_CharAttr_VType) != 'C')
            glvrvr[var] = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, 'C');
        var++;
    }
    model->update();
    for (auto& var2glvrvar: glvrvr) {
        model->addConstr( var2glvrvar.second <= U*(*var2glvrvar.first) );
        model->addConstr( var2glvrvar.second >= L*(*var2glvrvar.first) );
        model->addConstr( *u - U*(1 - (*var2glvrvar.first)) <= var2glvrvar.second );
        model->addConstr( *u - L*(1 - (*var2glvrvar.first)) >= var2glvrvar.second );
    }
}

void Gloverizer::linearize_qconstr(GRBQConstr* qconstr) {
    GRBQConstr* qc { qconstr };
    GRBQuadExpr qexpr { model->getQCRow(*qc) };
    GRBLinExpr new_lhs { qexpr.getLinExpr() };
    GRBVar* binvar { nullptr };
    for (unsigned int i = 0; i < qexpr.size(); i++) {
        GRBVar var1 { qexpr.getVar1(i) };
        GRBVar var2 { qexpr.getVar2(i) };
        if (var1.get(GRB_CharAttr_VType) != 'C')
            binvar = &var1;
        else
            binvar = &var2;
        for (auto& var2glvrvr: glvrvr)
            if (binvar->sameAs(*var2glvrvr.first))
                new_lhs += qexpr.getCoeff(i) * var2glvrvr.second;
    }
    model->addConstr(new_lhs,
                     qc->get(GRB_CharAttr_QCSense),
                     qc->get(GRB_DoubleAttr_QCRHS));
    model->remove(*qc);
}

void Gloverizer::linearize_objective() {
    GRBQuadExpr oldobj { model->getObjective() };
    GRBLinExpr newobj { oldobj.getLinExpr() };
    GRBVar* binvar { nullptr };
    for (unsigned int i = 0; i < oldobj.size(); i++) {
        GRBVar var1 { oldobj.getVar1(i) };
        GRBVar var2 { oldobj.getVar2(i) };
        if (var1.get(GRB_CharAttr_VType) != 'C')
            binvar = &var1;
        else
            binvar = &var2;
        for (auto& var2glvrvr: glvrvr)
            if (binvar->sameAs(*var2glvrvr.first))
                newobj += oldobj.getCoeff(i) * var2glvrvr.second;
    }
    model->setObjective(newobj, model->get(GRB_IntAttr_ModelSense));
}


} // namespace grbfrc
