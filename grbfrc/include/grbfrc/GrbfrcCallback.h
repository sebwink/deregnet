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

#ifndef GRBFRC_CALLBACK_H
#define GRBFRC_CALLBACK_H

#include <map>
#include <vector>
#include <memory>

#include <gurobi_c++.h>

namespace grbfrc {

template <typename T>
class Callback {

public:

    GRBCallback* yield(std::vector<GRBVar*>& tx,
                       std::vector<GRBVar*>* vars,
                       std::vector<std::map<T,GRBVar>>* vargrps = nullptr) {
        // check whether referenced variables are discrete
        if (vargrps) {
            transformed_vargrps = new std::vector<std::map<T, GRBVar>>();
            for (std::map<T, GRBVar>& vargrp : *(vargrps)) {
                std::map<T, GRBVar> transformed_vargrp;
                for (auto& e: vargrp)
                    transformed_vargrp[e.first] = *tx[getIndex(e.second, vars)];
                transformed_vargrps->push_back(transformed_vargrp);
            }
            this->register_vargrps();
        }
        return init_callback();
    }

    GRBCallback* yield(std::vector<std::map<T,GRBVar>>* vargrps = nullptr) {
        if (vargrps) {
            transformed_vargrps = new std::vector<std::map<T, GRBVar>>();
            *transformed_vargrps = *vargrps;
            this->register_vargrps();
        }
        return init_callback();
    }

    virtual void register_vargrps() {

    }

    GRBCallback* init_callback() {
        return this->cbp;
    }

protected:

    GRBCallback* cbp { nullptr };
    std::vector<std::map<T,GRBVar>>* transformed_vargrps { nullptr };

private:

    int getIndex(GRBVar& var, std::vector<GRBVar*>* vars) {
        for (unsigned int i = 0; i < vars->size(); i++)
            if (var.sameAs(*((*vars)[i]))) return i;
        return -1;
    }

};

template <typename T>
class GrbfrcCallback {

private:

    std::vector<std::map<T,GRBVar>>* vargrps;
    Callback<T>* cb;

public:

    GrbfrcCallback(Callback<T>* xcb, std::vector<std::map<T,GRBVar>>* xvargrps = nullptr)

        : vargrps { xvargrps }, cb { xcb } { }

    GRBCallback* yield(std::vector<GRBVar*>& tx,
                       std::vector<GRBVar*>* vars) {
        return cb->yield(tx, vars, vargrps);
    }

    GRBCallback* yield() {
        return cb->yield(vargrps);
    }

};

}    //     namespace deregnet


#endif
