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

#ifndef YGGY_H
#define YGGY_H

#include <string>
#include <map>

#include <gurobi_c++.h>
#include <grbfrc/FMILP.h>


/*
 *  max (cc*xc + cd*xd + d)/(ec*xc + ed*xd + f)
 *
 *  s.t.
 *
 *      Ac*xc + Ad*xd = b
 *      xc in R
 *      xd in Z
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

namespace grbfrc
{

class YGGY
 {
   private:

      FMILP* fmip;
      FMILPSol solution;
      // GrbFrcCallback* callback ...
      bool invalid;
      bool transformed;
      GRBModel transformation;
      GRBVar u;                    // u = 1 / (ec*xc + ed*xd + f)
      std::vector<GRBVar*> tx;     // tx ~ (u * xc, xd)

   public:


      YGGY(FMILP* fmipPtr);
      void printInvalidity();
      int transform();
      GRBModel getTransform();
      void solveTransform();
      void solveTransform(GRBCallback& callback);
      void run(int time_limit);
      void run(GRBCallback& callback, int time_limit);
      void writeSolution();
      FMILPSol getSolution();


   private:

      void backTransformSolution();
      int getIndex(GRBVar& var);
      void define_tx();
      void define_objective(GRBLinExpr& objNumerator, int objSense);
      void define_constraints(GRBLinExpr& objDenominator);
      bool getUmax(double& umax);

 };

}  // namespace grbfrc

#endif  // YGGY_H
