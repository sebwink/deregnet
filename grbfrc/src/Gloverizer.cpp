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
#ifdef GRBFRC_OMP
    omp_lock_t model_lock;
    omp_init_lock(&model_lock);
  #pragma omp parallel for
#endif
    for (int q = 0; q < model->get(GRB_IntAttr_NumQConstrs); q++) {
#ifndef GRBFRC_OMP
        linearize_qconstr(qconstr);
        qconstr++;
#else
        std::cout << q << std::endl;
        linearize_qconstr(qconstr, q, &model_lock);
#endif
    }
#ifdef GRBFRC_OMP
    omp_destroy_lock(&model_lock);
    for (int q = 0; q < model->get(GRB_IntAttr_NumQConstrs); q++) {
        model->remove(*qconstr);
        qconstr++;
    }
#endif

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

#ifndef GRBFRC_OMP
void Gloverizer::linearize_qconstr(GRBQConstr* qconstr) {
#else
void Gloverizer::linearize_qconstr(GRBQConstr* qconstr, int q, omp_lock_t* model_lock) {
#endif
#ifndef GRBFRC_OMP
    GRBQConstr* qc { qconstr };
#else
    GRBQConstr* qc { qconstr + q };
#endif
#ifdef GRBFRC_OMP
    omp_set_lock(model_lock);
#endif
    GRBQuadExpr qexpr { model->getQCRow(*qc) };
#ifdef GRBFRC_OMP
    omp_unset_lock(model_lock);
#endif
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
#ifdef GBRFRC_OMP
    omp_set_lock(model_lock);
#endif
    model->addConstr(new_lhs,
                     qc->get(GRB_CharAttr_QCSense),
                     qc->get(GRB_DoubleAttr_QCRHS));
#ifndef GRBFRC_OMP
    model->remove(*qc);
#endif
#ifdef GRBFRC_OMP
    omp_unset_lock(model_lock);
#endif
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
