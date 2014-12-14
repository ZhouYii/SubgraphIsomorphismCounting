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

#define PRINT_ON

// Governs whether subgraph mappings are printed out to file
//#define PRINT_SUBGRAPH_MATCHED
//#define SUBGRAPH_OUTPUT_FILE "matched_ouput.txt"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
#include "graph.cpp"
#include "turboiso.h"
#include "combinators.cpp"

/*
 * TYPEDEF
 */
typedef std::unordered_map<UINT, QueryNode*> query_node_map;

// For aggregating nodes in NEC Tree by label
typedef std::unordered_map<UINT, vector<QueryNode*>*> NECBucketMap;

typedef vector<vector<QueryNode*>*> VectorOfNeighborNodeVectors;

// NECNode->Label->Vector<Candidates>
typedef unordered_map<NECNode*, unordered_map<VERTEX, vector<VERTEX>*>* > CandidateRegions;

typedef unordered_map<VERTEX, vector<VERTEX>*> SubRegions;

typedef priority_queue<pair<NECNode*, UINT>*, vector<pair<NECNode*, UINT>*>, CompareNEC> NEC_PQ;
#endif
