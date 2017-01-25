#include <gurobi_c++.h>

#include <grbfrc/ObjVarTransform.h>
#include <grbfrc/Gloverizer.h>

namespace grbfrc {

ObjVarTransform::ObjVarTransform(FMILP* FmipPtr, double* objub, double* objlb)

                                 : fmip { FmipPtr },
                                   invalid { false },
                                   transformed { false },
                                   transformation { GRBModel(*(fmipPtr->baseModel)) },
                                   umax { objub },
                                   umin { objlb } {
    fmip->update();
}

void ObjVarTransform::printInvalidity() {
   // if (no umax available) "no umax"
   // if (no umin available) "no umin"
   // else if (continuous vars in denominator) "continuous variables in denominator"
}

int ObjVarTransform::transform() {
    // objective variable u = (ax + by + c) / (dx + ey + f)
    if (!getUmax()) return 127;
    if (!getUmin()) return 128;
    u = transformation.addVar(umin, umax, 0.0, GRB_CONTINUOUS);
    transformation.update();
    defineUConstr();
    transformation.update();
    Gloverizer gloverizer(&transformation, &u, umax);
    gloverizer.gloverize();
    return 0;
}



}         // namespace grbfrc
