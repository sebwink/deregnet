#ifndef USINGLEMON_H
#define USINGLEMON_H

#include <lemon/smart_graph.h>
#include <lemon/adaptors.h>

namespace deregnet {

#define INVALID lemon::INVALID

using Graph = lemon::SmartDigraph;

using Node = Graph::Node;

using NodeIt = Graph::NodeIt;

using OutArcIt = Graph::OutArcIt;
using InArcIt = Graph::InArcIt;

template<typename ValueType>
using NodeMap = Graph::NodeMap<ValueType>;

using NodeFilter = NodeMap<bool>;

using InducedSubgraph = lemon::FilterNodes<Graph>;


}  // namespace deregnet

#endif  // USINGLEMON_H
