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

#ifndef OBJVARTRANSFORM_H
#define OBJVARTRANSFORM_H

#include <vector>

#include <gurobi_c++.h>

#include <grbfrc/common.h>
#include <grbfrc/Algorithm.h>
#include <grbfrc/GrbfrcCallback.h>

/*
 *  max (cc*xc + cd*xd + d)/(ec*xc + ed*xd + f)
 *
 *  s.t.
 *
 *      Ac*xc + Ad*xd = b
 *      xc in R
 *      xd in Z
 */

// if ex != 0: Objective variable transform not applicable

/*
 *  u = (cc*xc + cd*xd + d)/(ed*xd + f)
 *
 */

/*
 *   max u
 *
 *   s.t.
 *
 *       Ac*xc + Ad*xd = b
 *       ed*(u*xd) + f*u = cc*xc + cd*xd + d
 *       xc, u in R
 *       xd in Z
 *
 */

// Gloverize( u * xd ) --> w = u * xd

/*
 *   max u
 *
 *   s.t.
 *
 *       Ac*xc + Ad*xd = b
 *       ed*w + f*u = cc*xc + cd*xd + d
 *       GloverConstraints(w, xd, u)
 *       xc, u in R
 *       xd in Z
 *
 */

namespace grbfrc {

class ObjVarTransform : public _Algorithm {

    private:

        FMILPSol solution;
        bool invalid;
        bool transformed;
        GRBModel transformation;
        GRBVar u;
        std::vector<GRBVar*> tx;
        double* umax;
        double* umin;

    public:

        ObjVarTransform(GRBEnv* env, double* objub = nullptr, double* objlb = nullptr);
        void printInvalidity();

        template <typename T>
        int transform(GrbfrcCallback<T>* cb) {

            copy_model();

            u = transformation.addVar(*umin, *umax, 0.0, GRB_CONTINUOUS);
            transformation.setObjective(GRBLinExpr(u), base_model->get(GRB_IntAttr_ModelSense));
            transformation.update();
            add_u_constraint();

            if (cb) {
                callback = cb->yield(tx, vars);
                transformation.setCallback(callback);
            }

            transformation.update();
            return 0;
        }

        GRBModel getTransform();
        void solveTransform();

        template <typename T>
        void run(GrbfrcCallback<T>* cb) {
            if (invalid)
                printInvalidity();
            else {
                std::cout << "\n=========== solving FMIP via Objective Variable Transform transform ===========\n\n";
                if (!transformed) transform(cb);
                solveTransform();
                backTransformSolution();
            }
        }

        void writeSolution(FMILPSol** xsolution);
        // FMILPSol getSolution();

    private:

        void copy_model();

        void backTransformSolution();
        int getIndex(GRBVar& var);
        void add_u_constraint();


};      // class ObjVarTransform

}       // namespace grbfrc

#endif  // OBJVARTRANSFORM_H
