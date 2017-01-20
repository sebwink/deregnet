// --------------------------------------------------------------------------
//                   deregnet -- find deregulated pathways
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

#include <set>

#include <deregnet/usinglemon.h>
#include <deregnet/StartHeuristic.h>

namespace deregnet {

Node* StartHeuristic::get_best_root() {
    double best;
    Node* best_node { nullptr };
    if (receptors)
        for (auto v : *receptors)
            update_best_node(&best_node, &v, &best);
    else
        for (NodeIt v(*graph); v != INVALID; ++v)
            update_best_node(&best_node, &v, &best);
    return best_node;
}


bool StartHeuristic::search_further() {
    size--;
    if (size > 0 && found_node)
        return true;
    else if (size > 0 && !found_node)
        return false;
    else if (size == 0) {
        success = true;
        return false;
    }
    else
        return false;
}

bool StartHeuristic::feasible_node(Node* node) {
    if (exclude) {
        if (start_solution->find(*node) == start_solution->end() && exclude->find(*node) == exclude->end())
            return true;
        }
    else
        if (start_solution->find(*node) == start_solution->end())
            return true;
    return false;
}

}    //    namespace deregnet
