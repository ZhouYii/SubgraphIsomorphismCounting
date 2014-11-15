#ifndef DEFINES_H
#define DEFINES_H

/*
 *  PARAMETERS
 * */
#define STARTING_VERTEX_CANDIDATE_SIZE 5

/*
 *  CONSTANT
 * */
#define UINT unsigned int
#define VERTEX unsigned int
#define LABEL unsigned int
#define ULONG unsigned long

#include <unordered_map>
#include <vector>
#include <algorithm>
#include "graph.cpp"

/*
 * TYPEDEF
 */
typedef std::unordered_map<UINT, QueryNode*> query_node_map;

// For aggregating nodes in NEC Tree by label
typedef std::unordered_map<UINT, vector<QueryNode*>*> NECBucketMap;

typedef vector<vector<QueryNode*>*> VectorOfNeighborNodeVectors;

#endif
