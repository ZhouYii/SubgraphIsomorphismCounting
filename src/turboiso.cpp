#ifndef TURBOISO_CPP
#define TURBOISO_CPP

#include <iostream>
#include <queue>        // std::priority_queue
#include <boost/foreach.hpp>
#include <unordered_set>

#include "parsefs.cpp"
#include "graph.cpp"
#include "defines.h"
#include "turboiso.h"

using namespace std;

QueryNode* GetStartingQueryVertex(DataGraph& dg, query_node_map* qg) {
    QueryNode* tmp;
    RankingPair* rp;
    priority_queue<RankingPair*, vector<RankingPair*>, CompareRankingPair> pq;

    /*
     * TODO : When the filters get implemented, refactor this
     */
    for(int index = 0; index < qg->size(); index += 1) {
        tmp = qg->at(index);
        double score_denominator = tmp->parents->size() + tmp->children->size();
        double score_numerator = dg.vfs->NumVertices(tmp->label);

        rp = new RankingPair;
        rp->v = tmp;
        rp->score = score_numerator / score_denominator;

        // Adjust PQ, keeping within fixed candidate size
        if(pq.size() < STARTING_VERTEX_CANDIDATE_SIZE) {
            pq.push(rp);
        } else {
            if(pq.top()->score > rp->score) {
                // TODO : Fix memory leak. Pop does not free the allocated
                // RankingPair
                delete pq.top();
                pq.pop();
                pq.push(rp);
            } else {
                delete rp;
            }
        }
    }
    tmp = pq.top()->v;
    while(pq.size() > 0) {
        rp = pq.top();
        pq.pop();
        delete rp;
    }

    /*
     * TODO : Two filters and CR estimation later
     */
    return tmp;
}

void GetSortedNeighborLabels(vector<QueryNode*>* nodes, vector<LABEL>* labels) {
    QueryNode* adj_node;
    for(int i = 0; i < nodes->size(); i += 1) {
        adj_node = nodes->at(i);
        labels->push_back(adj_node->label);
    }
    // Sort Ascending. Order doesn't actually matter.
    sort(labels->begin(), labels->end());
}

vector<NECNode*>* MatchNECNodes(NECBucketMap* node_buckets) {
    vector<NECNode*>* nec_nodes;
    vector<QueryNode*>* query_nodes;
    vector<LABEL> nec_child_labels;
    vector<LABEL> nec_parent_labels;
    vector<LABEL> cmp_child_labels;
    vector<LABEL> cmp_parent_labels;
    QueryNode* origin_node;
    QueryNode* cmp_node;
    NECNode* nec_node;

    nec_nodes = new vector<NECNode*>;
    for(NECBucketMap::iterator it = node_buckets->begin(); 
        it != node_buckets->end(); ++it) {
        query_nodes = it->second;

        while(query_nodes->size() > 0) {
            origin_node = query_nodes->back();
            query_nodes->pop_back();
            nec_node = NECNodeInit(); // Allocates memory for member variables as well
            nec_node->members->push_back(origin_node);
            nec_nodes->push_back(nec_node);

            // Get adjacent information
            // This is the slow path. Take this after the hash path :
            // TODO : two separate hashes for parent and child labels
            GetSortedNeighborLabels(origin_node->children, &nec_child_labels);
            GetSortedNeighborLabels(origin_node->parents, &nec_parent_labels);

            // Iterate over all remaining nodes in vector (in reverse) to find
            for(int node_index = query_nodes->size() - 1;
                node_index >= 0; node_index -= 1) {
                cmp_node = query_nodes->at(node_index);
                GetSortedNeighborLabels(cmp_node->children, &cmp_child_labels);
                GetSortedNeighborLabels(cmp_node->parents, &cmp_parent_labels);
                if (cmp_child_labels == nec_child_labels &&
                    cmp_parent_labels == nec_parent_labels) {
                    query_nodes->erase(query_nodes->begin() + node_index);
                    nec_node->members->push_back(cmp_node);
                }
                cmp_parent_labels.clear();
                cmp_child_labels.clear();
            }
            // same adjacent signatures
            nec_child_labels.clear();
            nec_parent_labels.clear();
        }
    }

    return nec_nodes;
}

vector<NECNode*>* GenerateNECNodes(unordered_set<QueryNode*>* pnodes,
                                   unordered_set<QueryNode*>* seen)
{
    vector<QueryNode*>* bucket_vector;
    NECBucketMap query_node_buckets;
    QueryNode* curr_node;
    
    unordered_set<QueryNode*>::iterator it = pnodes->begin();
    for(; it != pnodes->end(); ++it) {
        curr_node = *(it);

        // Check if node has been marked.
        if(seen->count(curr_node) == 1) {
            continue;
        }
        // Mark node
        seen->insert(curr_node);
        if (query_node_buckets.count(curr_node->label) == 0) {
            query_node_buckets.insert(
                make_pair<LABEL,vector<QueryNode*>*>(curr_node->label,new vector<QueryNode*>));
        }
        bucket_vector = query_node_buckets[curr_node->label];
        bucket_vector->push_back(curr_node);
    }

    return MatchNECNodes(&query_node_buckets);
}

unordered_set<QueryNode*>* GetAdjacentChildQueryNodes(NECNode* n) {
    QueryNode* curr;
    unordered_set<QueryNode*>* s = new unordered_set<QueryNode*>;
    for (int member = 0; member < n->members->size(); member += 1) {
        curr = n->members->at(member);
        for(int child = 0; child < curr->children->size(); child += 1) {
            s->insert(curr->children->at(child));
        }
    }
    return s;
}

unordered_set<QueryNode*>* GetAdjacentParentQueryNodes(NECNode* n) {
    QueryNode* curr;
    unordered_set<QueryNode*>* s = new unordered_set<QueryNode*>;
    for (int member = 0; member < n->members->size(); member += 1) {
        curr = n->members->at(member);
        for(int parent = 0; parent < curr->parents->size(); parent += 1) {
            s->insert(curr->parents->at(parent));
        }
    }
    return s;
}

NECNode* RewriteToNECTree(query_node_map* qg, QueryNode* start, 
                          unordered_set<NECNode*>* leaf_nec_nodes) {
    NECNode* curr;
    NECNode* root;
    vector<NECNode*>* curr_level_nec = new vector<NECNode*>;
    vector<NECNode*>* next_level_nec = new vector<NECNode*>;
    vector<NECNode*>* tmp;
    vector<QueryNode*> candidate_set;
    unordered_set<QueryNode*> visited;
    unordered_set<QueryNode*>* adj_query_nodes;

    curr = NECNodeInit();
    root = curr;
    curr->members->push_back(start);
    visited.insert(start);
    next_level_nec->push_back(curr);
    do {
        // Swap curr and next literation
        tmp = curr_level_nec;
        curr_level_nec = next_level_nec;
        next_level_nec = tmp;
        tmp->clear();

        for (int nec_idx = 0; nec_idx < curr_level_nec->size(); nec_idx += 1) {
            curr = curr_level_nec->at(nec_idx);

            // Set of adj child QueryNodes
            adj_query_nodes = GetAdjacentChildQueryNodes(curr);
            vector<NECNode*>* children = GenerateNECNodes(adj_query_nodes, &visited);
            curr->children = children;
            for (int i = 0; i < children->size(); i += 1)
                next_level_nec->push_back(children->at(i));
            delete adj_query_nodes;

            // Set of adj parent QueryNodes
            adj_query_nodes = GetAdjacentParentQueryNodes(curr);
            vector<NECNode*>* parents = GenerateNECNodes(adj_query_nodes, &visited);
            curr->parents = parents;
            for (int i = 0; i < parents->size(); i += 1)
                next_level_nec->push_back(parents->at(i));
            delete adj_query_nodes;

            if (curr->children->size() == 0) {
                leaf_nec_nodes->insert(curr);
            }

            cout << "NEC Node ID :" << curr << " Children:";
            for (int i = 0; i < curr->children->size(); i += 1) {
                cout << curr->children->at(i) << " | ";
            }
            cout << " Parents:";
            for (int i = 0; i < curr->parents->size(); i += 1) {
                cout << curr->parents->at(i) << " | ";
            }
            cout << endl;
        }
    } while (curr_level_nec->size() > 0);

    return root;
}

bool DegreeFilter(DataGraph& dg, VERTEX data_node_id, NECNode* nec_root) {
    QueryNode* nec_node = nec_root->members->at(0);
    return nec_node->children->size() == dg.efs->GetOutDegree(data_node_id);
}

/* TODO : Sort adjacency list of query nodes by ascending label so we can do
 * this quickly
bool NeighborLabelFilter(DataGraph& dg, UINT data_node_id, NECNode* nec_root) {
    QueryNode* nec_node = nec_root->members->at(0);

}
*/

/*
bool ExploreCandidateRegions(DataGraph& dg, NECNode* nec_root,
                             vector<VERTEX>* data_node_candidates,
                             unordered_set<VERTEX>* marked_vertex, CR_ROOT)
{
    vector<VERTEX>::iterator it;
    VERTEX data_node_id;
    bool matched;

    it = data_node_candidates->begin();
    for(; it != data_node_candidates->end(); ++it) {
        data_node_id = *(it);
        if (marked_vertex->count(data_node_id) > 0 ||
            DegreeFilter(dg, data_node_id, nec_root) == true ||
            NeighborLabelFilter(dg, data_node_id, nec_root))
            continue;

        marked_vertex->insert(data_node_id);
        matched = true;

        // Need to sort children

    }

    if (CR_ROOT.size() < nec_root->members->size()) {
        ClearAllCR();
        return false;
    }
    return true;
}
*/

void TurboIso(DataGraph& dg, query_node_map* qg) {
    QueryNode* start_vertex;
    VERTEX* data_node_iter;
    LABEL root_label;
    VERTEX root_node_id;
    unordered_set<VERTEX>* visited_vertices = new unordered_set<VERTEX>;
    unordered_set<NECNode*>* leaf_nec_nodes = new unordered_set<NECNode*>;
    vector<VERTEX> candidate_data_vertices;
    NECNode* nec_root;

    start_vertex = GetStartingQueryVertex(dg, qg);
    root_label = start_vertex->label;
    nec_root = RewriteToNECTree(qg, start_vertex, leaf_nec_nodes);
    /*

    data_node_iter = dg.vfs->GetVertexIterator(root_label);
    for(int data_node_index = 0; 
        data_node_index < dg.vfs->NumVertices(root_label); 
        data_node_index += 1)
    {
        root_node_id = *(data_node_iter + data_node_index);
        candidate_data_vertices.clear();
        candidate_data_vertices.push_back(root_node_id);
        if(ExploreCandidateRegions(dg, nec_root, &candidate_data_vertices,
                                   visited_vertices, CR_ROOT) == false) 
        {
            continue;
        }
    }
    */
}

int main() {
    DataGraph dg;
    query_node_map* qg;

    VertexFs vfs("test.txt");
    EdgeFs efs("test_efs.txt", &vfs);
    dg.vfs = &vfs;
    dg.efs = &efs;

    qg = ReadQueryGraphFromFile("query_graph.txt");

    TurboIso(dg, qg);
}
#endif
