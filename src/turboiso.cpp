#ifndef TURBOISO_CPP
#define TURBOISO_CPP

#include <iostream>
#include <queue>        // std::priority_queue
#include <ctime>
#include <boost/foreach.hpp>
#include <unordered_set>

#include "parsefs.cpp"
#include "graph.cpp"
#include "defines.h"
#include "turboiso.h"

using namespace std;
bool ExploreCandidateRegions(DataGraph& dg, QueryGraph* qg,
                             NECNode* nec_root,
                             vector<VERTEX>* data_node_candidates,
                             unordered_set<VERTEX>* marked_vertex,
                             CandidateRegions* update_cr,
                             CandidateRegions* secondary_cr,
                             CandidateRegions* child_cr,
                             VERTEX parent,
                             NECNode* parent_nec,
                             unordered_set<NECNode*>* visited_nec);

pair<NECNode*, VERTEX>* MakePair(NECNode* nec, VERTEX prior) {
    pair<NECNode*, VERTEX>* p = new pair<NECNode*, VERTEX>;
    p->first = nec;
    p->second = prior;
    return p;
}

QueryNode* GetStartingQueryVertex(DataGraph& dg, QueryGraph* query_graph) {
    QueryNode* tmp;
    RankingPair* rp;
    priority_queue<RankingPair*, vector<RankingPair*>, CompareRankingPair> pq;
    query_node_map* qg = query_graph->query_nodes;

    // Any initialized query node should be it's own NEC. Minimizes CR size
    // since the NEC contains only one query node. 
    if (query_graph->initialized_query_nodes->size() > 0)
        return query_graph->initialized_query_nodes->begin()->first;

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
    return tmp;
}

/*
void GetSortedNeighbors(vector<QueryNode*>* nodes, vector<QueryNode*>* els) {
    for(int i = 0; i < nodes->size(); i += 1)
        labels->push_back(nodes->at(i));

    // Sort Ascending. Order doesn't actually matter.
    sort(labels->begin(), labels->end());
}
*/

void MatchNECNodes(NECBucketMap* node_buckets, vector<NECNode*>* buffer) {
    vector<QueryNode*>* query_nodes;
    vector<LABEL> nec_child_labels;
    vector<LABEL> nec_parent_labels;
    vector<LABEL> cmp_child_labels;
    vector<LABEL> cmp_parent_labels;
    QueryNode* query_node;
    QueryNode* origin_node;
    QueryNode* cmp_node;
    NECNode* nec_node;

    for(NECBucketMap::iterator it = node_buckets->begin(); 
        it != node_buckets->end(); ++it) {
        query_nodes = it->second;

        for (int index = 0; index < query_nodes->size(); index += 1) {
            query_node = query_nodes->at(index);
            sort(query_node->children->begin(), query_node->children->end());
            sort(query_node->parents->begin(), query_node->parents->end());
        }

        while(query_nodes->size() > 0) {
            origin_node = query_nodes->back();
            query_nodes->pop_back();
            nec_node = NECNodeInit(); // Allocates memory for member variables as well
            nec_node->members->push_back(origin_node);
            buffer->push_back(nec_node);

            for(int node_index = query_nodes->size() - 1;
                node_index >= 0; node_index -= 1) 
            {
                cmp_node = query_nodes->at(node_index);
                // Not sorting labels but memory addresses. This ensures same
                // adjacent query nodes in each NEC
                if (*(origin_node->children) == *(cmp_node->children) &&
                    *(origin_node->parents) == *(cmp_node->parents)) {
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
}

void GenerateNECNodes(unordered_set<QueryNode*>* pnodes,
                                   unordered_set<QueryNode*>* seen,
                                   vector<NECNode*>* buffer)
{
    vector<QueryNode*>* bucket_vector;
    NECBucketMap query_node_buckets;
    QueryNode* curr_node;
    
    unordered_set<QueryNode*>::iterator it = pnodes->begin();
    // Group Query Nodes by NEC
    for(; it != pnodes->end(); ++it) {
        curr_node = *(it);

        // Check if node has been marked.
        if(seen->count(curr_node) > 0) {
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

    MatchNECNodes(&query_node_buckets, buffer);
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

NECNode* RewriteToNECTree(QueryGraph* query_graph, 
                          QueryNode* start, 
                          vector<NECNode*>& initialized_nec) {
    NECNode* curr;
    NECNode* root;
    vector<NECNode*>* curr_level_nec = new vector<NECNode*>;
    vector<NECNode*>* next_level_nec = new vector<NECNode*>;
    vector<NECNode*>* tmp;
    vector<QueryNode*> candidate_set;
    unordered_set<QueryNode*> visited;
    unordered_set<NECNode*> visited_nec;
    unordered_set<QueryNode*>* adj_query_nodes;
    query_node_map* qg = query_graph->query_nodes;
    unordered_map<QueryNode*, VERTEX>::iterator initialized_it;

    if (query_graph->initialized_query_nodes->size() == 0) {
        // No initialized vertex Case
        root = NECNodeInit();
        root->members->push_back(start);
        visited.insert(start);
        visited_nec.insert(root);
        next_level_nec->push_back(root);
    } else {
        // Use a supernode for NEC when there are initialized query vertices
        root = NECNodeInit();
        root->members->clear();
        visited_nec.insert(root);
        initialized_it = query_graph->initialized_query_nodes->begin();
        while (initialized_it != query_graph->initialized_query_nodes->end()) {
            visited.insert(initialized_it->first);
            curr = NECNodeInit();
            visited_nec.insert(curr);
            curr->members->push_back(initialized_it->first);
            // curr->parents->push_back(root);  // Do not form edge to supernode
            next_level_nec->push_back(curr);
            root->children->push_back(curr);
            initialized_nec.push_back(curr);
            initialized_it++;
        }
    }

    do {
        // Swap curr and next iteration
        tmp = curr_level_nec;
        curr_level_nec = next_level_nec;
        next_level_nec = tmp;
        next_level_nec->clear();

        for (int nec_idx = 0; nec_idx < curr_level_nec->size(); nec_idx += 1) {
            curr = curr_level_nec->at(nec_idx);
            visited_nec.insert(curr);

            // Set of adj child QueryNodes
            adj_query_nodes = GetAdjacentChildQueryNodes(curr);
            GenerateNECNodes(adj_query_nodes, &visited, curr->children);
            for (int i = 0; i < curr->children->size(); i += 1) {
                if (visited_nec.count(curr->children->at(i)) == 0) {
                    curr->children->at(i)->parents->push_back(curr);
                    next_level_nec->push_back(curr->children->at(i));
                }
            }
            delete adj_query_nodes;

            // Set of adj parent QueryNodes
            adj_query_nodes = GetAdjacentParentQueryNodes(curr);
            GenerateNECNodes(adj_query_nodes, &visited, curr->parents);
            for (int i = 0; i < curr->parents->size(); i += 1) {
                if (visited_nec.count(curr->parents->at(i)) == 0) {
                    curr->parents->at(i)->children->push_back(curr);
                    next_level_nec->push_back(curr->parents->at(i));
                }
            }
            delete adj_query_nodes;
        }
    } while (curr_level_nec->size() > 0);

    delete curr_level_nec;
    delete next_level_nec;

    return root;
}

bool DegreeFilter(DataGraph& dg, VERTEX data_node_id, NECNode* nec_root) {
    QueryNode* nec_node = nec_root->members->at(0);
    bool child_mismatch = nec_node->children->size() > dg.efs->GetOutDegree(data_node_id);
    bool parent_mismatch =  nec_node->parents->size() > dg.efs->GetInDegree(data_node_id);
    return child_mismatch || parent_mismatch;
}

// Since CandidateRegions(cr, nec, data_vertex) is deterministic operator, we
// don't need to ever calculate candidate regions for a nec/data_vertex pair
// twice.
bool CandidateRegionBeenVisited(CandidateRegions* cr, 
                                NECNode* nec, 
                                VERTEX prior,
                                unordered_set<NECNode*>* visited_nec) 
{
    if (visited_nec->count(nec) > 0) {
        return true;
    }

    unordered_map<VERTEX, vector<VERTEX>*>* prior_vertex_map;
    if (cr->count(nec) != 0) {
        prior_vertex_map = cr->find(nec)->second;
        if (prior_vertex_map->count(prior) != 0) {
            return true;
        }
    }
    return false;
}

/* TODO : Sort adjacency list of query nodes by ascending label so we can do
 * this quickly
bool NeighborLabelFilter(DataGraph& dg, UINT data_node_id, NECNode* nec_root) {
    QueryNode* nec_node = nec_root->members->at(0);

}
*/

/* Rank children NEC by number of candidate data vertices to form evaluation
 * order
 */
NEC_PQ* GetChildNECEvalOrder(DataGraph& dg, 
                             NECNode* nec_root,
                             VERTEX v_index, 
                             CandidateRegions* cr,
                             unordered_set<NECNode*>* visited_nec)
{
    vector<NECNode*>* nec_nodes;
    vector<NECNode*>::iterator it;
    pair<NECNode*, UINT>* p;
    NECNode* curr;
    NEC_PQ* pq = new NEC_PQ;

    cout << "child nec to visit : ";

    nec_nodes = nec_root->children;
    for (it = nec_nodes->begin(); it != nec_nodes->end(); it++) {
        curr = *(it);

        if (curr->members->size() == 0)
            continue;

        if (CandidateRegionBeenVisited(cr, curr, v_index, visited_nec))
            continue;

        cout << curr << " ";
        p = new pair<NECNode*, UINT>;
        p->first = curr;
        p->second = dg.efs->GetOutDegree(v_index, curr->members->at(0)->label);
        pq->push(p);
    }
    cout << endl;

    return pq;
}

NEC_PQ* GetParentNECEvalOrder(DataGraph& dg, 
                              NECNode* nec_root,
                              VERTEX v_index, 
                              CandidateRegions* cr,
                              unordered_set<NECNode*>* visited_nec) 
{
    vector<NECNode*>* nec_nodes;
    vector<NECNode*>::iterator it;
    pair<NECNode*, UINT>* p;
    NECNode* curr;
    NEC_PQ* pq = new NEC_PQ;

    cout << "parent nec to visit : ";

    nec_nodes = nec_root->parents;
    for (it = nec_nodes->begin(); it != nec_nodes->end(); it++) {
        curr = *(it);

        if (curr->members->size() == 0)
            continue;

        if (CandidateRegionBeenVisited(cr, curr, v_index, visited_nec))
            continue;

        cout << curr << " ";

        p = new pair<NECNode*, UINT>;
        p->first = curr;
        p->second = dg.efs->GetInDegree(v_index, curr->members->at(0)->label);
        pq->push(p);
    }

    cout << endl;
    return pq;
}

void ClearCandidateRegions(CandidateRegions* cr, 
                           vector<pair<NECNode*,VERTEX>*>* new_cr_entries)
{
    pair<NECNode*, VERTEX>* pair;
    for (int index = 0; index < new_cr_entries->size(); index += 1) {
        pair = new_cr_entries->at(index);
        cr->at(pair->first)->at(pair->second)->clear();
        delete pair;
    }
    /*
    NECNode* clear_nec;

    for (int neighbor = 0; neighbor < visited_nec_neighbors->size(); neighbor += 1) 
    {
        clear_nec = visited_nec_neighbors->at(neighbor);
        if (cr->count(clear_nec) > 0 &&
            cr->find(clear_nec)->second->count(parent) > 0) {
            cr->find(clear_nec)->second->find(parent)->second->clear();
        }
    }
    */
}

void CandidateRegionAddNEC(CandidateRegions* cr, NECNode* nec) 
{
    if (cr->count(nec) == 0) {
        cr->insert(make_pair<NECNode*, unordered_map<VERTEX, vector<VERTEX>*>*>
                        (nec, new unordered_map<VERTEX, vector<VERTEX>*>));
    }
}


void CandidateRegionAddPrior(CandidateRegions* cr,
                             NECNode* nec,
                             VERTEX prior,
                             vector<pair<NECNode*, VERTEX>*>* new_cr_entries) 
{
    unordered_map<VERTEX, vector<VERTEX>*>* prior_vertex_map;
    CandidateRegionAddNEC(cr, nec);
    prior_vertex_map = cr->find(nec)->second;
    if (prior_vertex_map->count(prior) == 0) {
        prior_vertex_map->insert(make_pair<VERTEX, vector<VERTEX>*>
                                          (prior, new vector<VERTEX>));
        new_cr_entries->push_back(new pair<NECNode*, VERTEX>(nec, prior));
    }
}

/*  Function adds a new entry to update_cr: (nec_node, parent) -> vertex_id
 *  Also adds new entry to secondary_cr: (parent_nec, vertex_id) -> parent
 *  Newly created cr entries are updated in the new_update_cr and
 *  new_secondary_cr structures.
 *
 *  Parameters :
 *  update_cr, secondary_cr :
 *      two candidate regions corresponding to either parent_cr/child_cr or vice
 *      child_cr/parent_cr depending how the function is called.
 *  parent : 
 *      the prior data node for the inserted cr entry
 *  parent_nec : 
 *      NEC class of the parent data node
 *  nec_node : 
 *      NEC class of the inserted data node to the cr
 *  vertex_id :
 *      The inserted data not to the cr
 *  new_update_cr, new_seconday_cr :
 *      new_update_cr corresponds to the newly created entries in update_cr
 *      These data structures help selectively clear candidate regions without
 *      clearing valid neighboring candidate regions which were discovered
 *      elsewhere in the recursion stack.
 */
void CandidateRegionUpdate(CandidateRegions* update_cr,
                           CandidateRegions* secondary_cr,
                           VERTEX parent,
                           NECNode* parent_nec,
                           NECNode* nec_node,
                           VERTEX vertex_id,
                           vector<pair<NECNode*, VERTEX>*>* new_update_cr,
                           vector<pair<NECNode*, VERTEX>*>* new_secondary_cr) 
{
    CandidateRegionAddPrior(update_cr, nec_node, parent, new_update_cr);
    CandidateRegionAddPrior(secondary_cr, parent_nec, vertex_id, new_secondary_cr);
    update_cr->find(nec_node)->second->find(parent)->second->push_back(vertex_id);
    secondary_cr->find(parent_nec)->second->find(vertex_id)->second->push_back(parent);
}

// Unfortunately we cannot pass pointers to memberfunctions of
// edgeFs->getIn/Outvertex
bool ExploreParentSubRegions(DataGraph& dg,
                       QueryGraph* qg,
                       unordered_set<VERTEX>* marked_vertex,  
                       CandidateRegions* update_cr,
                       CandidateRegions* secondary_cr,
                       CandidateRegions* child_cr,
                       VERTEX parent,
                       VERTEX adj_vertex,
                       NECNode* adj_vertex_nec,
                       vector<NECNode*>* visited_nec_neighbors,
                       NEC_PQ* pq,
                       vector<pair<NECNode*,VERTEX>*>* new_child_cr,
                       vector<pair<NECNode*,VERTEX>*>* new_parent_cr,
                       unordered_set<NECNode*>* visited_nec)
{
    NECNode* adj_nec;
    vector<UINT>* candidate_vertices;
    LABEL candidate_label;
    bool matched = true;
    CandidateRegions* parent_cr = update_cr;
    if (child_cr == update_cr)
        parent_cr = secondary_cr;

    if (pq->size() > 0) {
        // Child NEC
        adj_nec = pq->top()->first;
        delete pq->top();
        pq->pop();

        // To decide whether a ChildCR should be explored, given prior knowledge
        // that the adj_vertex data vertex will be mapped to adj_vertex_nec, 
        //
        // Exploring a child NEC vertex means a prior vertex for a parent NEC
        // has been decided. Therefore the Candidate Region that will be updated
        // will be the ParentCR(child_nec, prior) since candidate data vertices
        // for the child_nec must have prior as a parent data vertex.
        if (!CandidateRegionBeenVisited(parent_cr, adj_nec, adj_vertex, visited_nec)) {
            candidate_label = adj_nec->members->at(0)->label;
            // TODO : this is a bug, since this is not able of duality between
            // parent and child.
            candidate_vertices = dg.efs->GetIndegreeVertices(adj_vertex, 
                                                            candidate_label);
#ifdef PRINT_ON
            cout << "candidates for " << adj_nec << ": ";
            for (int i = 0; i < candidate_vertices->size(); i += 1) {
                cout << candidate_vertices->at(i) << " - ";
            }
            cout << endl;
#endif
            // Construct candidate vertex using VFS
            if (!ExploreCandidateRegions(dg, qg, adj_nec, candidate_vertices,
                                        marked_vertex, update_cr, secondary_cr,  child_cr,
                                        adj_vertex, adj_vertex_nec, visited_nec))
            {
                ClearCandidateRegions(child_cr, new_child_cr);
                ClearCandidateRegions(parent_cr,new_parent_cr);
                matched = false;
            }
            visited_nec_neighbors->push_back(adj_nec);
            delete candidate_vertices;
        }
    }
    return matched;
}

// Attept to match adj_vertex to adj_vertex_nec
bool ExploreChildSubRegions(DataGraph& dg,
                       QueryGraph* qg,
                       unordered_set<VERTEX>* marked_vertex,  
                       CandidateRegions* update_cr,
                       CandidateRegions* secondary_cr,
                       CandidateRegions* child_cr,
                       VERTEX parent,
                       VERTEX adj_vertex,
                       NECNode* adj_vertex_nec,
                       vector<NECNode*>* visited_nec_neighbors,
                       NEC_PQ* pq,
                       vector<pair<NECNode*,VERTEX>*>* new_child_cr,
                       vector<pair<NECNode*,VERTEX>*>* new_parent_cr,
                       unordered_set<NECNode*>* visited_nec)
{
    NECNode* adj_nec;
    vector<UINT>* candidate_vertices;
    LABEL candidate_label;
    bool matched = true;
    CandidateRegions* parent_cr = update_cr;
    if (child_cr == update_cr)
        parent_cr = secondary_cr;

    if (pq->size() > 0) {
        // Child NEC
        adj_nec = pq->top()->first;
        delete pq->top();
        pq->pop();

        // To decide whether a ChildCR should be explored, given prior knowledge
        // that the adj_vertex data vertex will be mapped to adj_vertex_nec, 
        //
        // Exploring a child NEC vertex means a prior vertex for a parent NEC
        // has been decided. Therefore the Candidate Region that will be updated
        // will be the ParentCR(child_nec, prior) since candidate data vertices
        // for the child_nec must have prior as a parent data vertex.
        if (!CandidateRegionBeenVisited(parent_cr, adj_nec, adj_vertex, visited_nec)) {
            candidate_label = adj_nec->members->at(0)->label;
            // TODO : this is a bug, since this is not able of duality between
            // parent and child.
            candidate_vertices = dg.efs->GetOutdegreeVertices(adj_vertex, 
                                                            candidate_label);
#ifdef PRINT_ON
            cout << "candidates for " << adj_nec << ": ";
            for (int i = 0; i < candidate_vertices->size(); i += 1) {
                cout << candidate_vertices->at(i) << " - ";
            }
            cout << endl;
#endif
            // Construct candidate vertex using VFS
            if (!ExploreCandidateRegions(dg, qg, adj_nec, candidate_vertices,
                                        marked_vertex, update_cr, secondary_cr, child_cr,
                                        adj_vertex, adj_vertex_nec, visited_nec))
            {
                ClearCandidateRegions(child_cr, new_child_cr);
                ClearCandidateRegions(parent_cr,new_parent_cr);
                matched = false;
            }
            visited_nec_neighbors->push_back(adj_nec);
            delete candidate_vertices;
        }
    }
    return matched;
}

// Verify that the number of candidates to a candidate region is at least the
// size of the number of query vertices of the NEC
bool CandidateRegionVerifySize(CandidateRegions* cr, VERTEX prior, NECNode* nec) {
    bool ret;
    if (cr->count(nec) == 0)
        ret = false;
    else {
        ret = true;

        if (cr->find(nec)->second->count(prior) == 0) {
            ret = false;
        } else {
            if (cr->find(nec)->second->find(prior)->second->size() 
                            < nec->members->size() && ret == true) 
            {
                cr->find(nec)->second->find(prior)->second->clear();
                ret = false;
            } else {
                ret = true;
            }
        }
    }
    return ret;
}

// Updates Child Candidate Regions(nec_root, parent)
bool ExploreCandidateRegions(DataGraph& dg, 
                             QueryGraph* qg,
                             NECNode* nec_root,
                             vector<VERTEX>* data_node_candidates,
                             unordered_set<VERTEX>* marked_vertex,
                             CandidateRegions* update_cr,
                             CandidateRegions* secondary_cr,
                             CandidateRegions* child_cr,
                             VERTEX parent,
                             NECNode* parent_nec,
                             unordered_set<NECNode*>* visited_nec)
{
#ifdef PRINT_ON
    cout << "Explore CR: " << nec_root << endl;
    cout << "Data Candidates: ";
    for (int i = 0; i < data_node_candidates->size(); i += 1) {
        cout << data_node_candidates->at(i) << "-";
    }
    cout << endl;
#endif
    vector<VERTEX>::iterator it;
    NEC_PQ* child_nec_pq;
    NEC_PQ* parent_nec_pq;
    vector<NECNode*> visited_nec_neighbors;
    vector<pair<NECNode*, VERTEX>*> new_child_cr;
    vector<pair<NECNode*, VERTEX>*> new_parent_cr;
    VERTEX data_node_id;
    bool matched;
    bool ret;

    CandidateRegions* parent_cr = update_cr;
    if (child_cr == update_cr)
        parent_cr = secondary_cr;

    // Completed NEC won't be needed as a way of only visiting each NEC once -
    // each NEC may be visited many times in order to mark all pairs of
    // NEC,Prior. The marking process will be at the step of iterator. The
    // marked sets require parent/child separation. A data node prior as child
    // gives different candidate vertices compared to the same data node prior
    // as parent.
    cout << " ---- INSERT:" << nec_root << endl;
    visited_nec->insert(nec_root);
    it = data_node_candidates->begin();
    for(; it != data_node_candidates->end(); ++it) {
        // data_node_id is drawn from set of data node candidates which are
        // compatible with nec_root. Therefore, the nec_root is a compatible
        // NEC for the data_node_id.
        data_node_id = *(it);
        if (marked_vertex->count(data_node_id) > 0 ||
            DegreeFilter(dg, data_node_id, nec_root) == true)
            //NeighborLabelFilter(dg, data_node_id, nec_root))
            continue;

        marked_vertex->insert(data_node_id);
        matched = true;

        // Need to sort children by evaluation cost so that the data node can be
        // failed earlier for incompatiblity reasons.
        // We can also prune out children NEC + parent NEC by checking if
        // <datanode/parent> + NEC combination has already been visited. This
        // means we will not need to check during subregions exploration.
        child_nec_pq = GetChildNECEvalOrder(dg,
                                            nec_root,
                                            data_node_id,
                                            child_cr,
                                            visited_nec);
        parent_nec_pq = GetParentNECEvalOrder(dg,
                                              nec_root,
                                              data_node_id,
                                              parent_cr,
                                              visited_nec);

        // Recursively generate CR for child NEC, by interleaving parent and
        // child NEC nodes.
        while (child_nec_pq->size() > 0 || parent_nec_pq->size() > 0) {
            // Due to the children/parent pq generation step, this block of code
            // only executes of previously unvisited pairs of prior_data_vertex
            // and NEC node. When exploring child or parent NEC, need to mark
            // the data vertex and NEC pair as visited
            if (// Update child_cr
                !ExploreChildSubRegions(dg, qg, marked_vertex, parent_cr, child_cr,
                                  child_cr,
                                  parent, data_node_id, nec_root,
                                  &visited_nec_neighbors, child_nec_pq,
                                  &new_child_cr, &new_parent_cr, visited_nec) ||
                // Update parent_cr
                !ExploreParentSubRegions(dg, qg, marked_vertex, child_cr, parent_cr,
                                  child_cr,
                                  parent, data_node_id, nec_root,
                                  &visited_nec_neighbors, parent_nec_pq,
                                  &new_child_cr, &new_parent_cr, visited_nec)) {
                matched = false;
                break;
            }
        }

        marked_vertex->erase(data_node_id);
        if (!matched)
            continue;

        CandidateRegionUpdate(update_cr, secondary_cr, parent, parent_nec,
                              nec_root, data_node_id, &new_child_cr,
                              &new_parent_cr);
    }
    visited_nec->erase(nec_root);
    ret = CandidateRegionVerifySize(update_cr, parent, nec_root);
    return ret;
}

vector<NECNode*>* MatchingOrderGetUnvisited(NECNode* leaf, 
                                        unordered_set<NECNode*>* visited) {
    vector<NECNode*>* v = new vector<NECNode*>;
    NECNode* curr;
    for (int i = 0; i < leaf->children->size(); i += 1) {
        curr = leaf->children->at(i);
        if (visited->count(curr) == 0)
            v->push_back(curr);
    }
    for (int i = 0; i < leaf->parents->size(); i += 1) {
        curr = leaf->parents->at(i);
        if (visited->count(curr) == 0)
            v->push_back(curr);
    }
    return v;
}

void MatchingOrderVisitChild(NECNode* child, vector<NECNode*>* path, ULONG nontree_edges,
                             queue<pair<vector<NECNode*>*, ULONG>*>* bfs_queue,
                             unordered_set<NECNode*>* visited) {
    vector<NECNode*>* new_path;
    pair<vector<NECNode*>*, ULONG>* new_pair;

    if (visited->count(child) > 0) {
        // Already marked
        return;
    }

    visited->insert(child);
    new_path = new vector<NECNode*>(*path);

    new_path->push_back(child);
    new_pair = new pair<vector<NECNode*>*, ULONG>;
    new_pair->first = new_path;
    new_pair->second = nontree_edges;
    bfs_queue->push(new_pair);
}

ULONG multiplyRange(ULONG start, ULONG end) {
    if (start == end)
        return start;
    return multiplyRange(start, end/2) * multiplyRange(end/2 + 1, end);
}

ULONG xChoseY(ULONG x, ULONG y) {
    return multiplyRange(y+1, x);
}

ULONG MatchingOrderCRScore(CandidateRegions* cr,
                          ULONG nontree_edges,
                          NECNode* leaf) 
{
    unordered_map<VERTEX, vector<VERTEX>*>* cr_nodes;
    ULONG score = 0;
    unordered_map<VERTEX, vector<VERTEX>*>::iterator it;
    ULONG member_size = leaf->members->size();

    if (cr->count(leaf) > 0) {
        cr_nodes = cr->at(leaf);
        it = cr_nodes->begin();
        if (member_size > 1) {
            while (it != cr_nodes->end()) {
                score += xChoseY(it->second->size(), member_size);
                ++ it;
            }
        } else {
            while (it != cr_nodes->end()) {
                score += it->second->size();
                ++it;
            }
            score = score / (nontree_edges + 1);
        }
    }
    return score;
}

ULONG MatchingOrderCalcScore(vector<NECNode*>* path,
                             ULONG nontree_edges,
                             CandidateRegions* child_cr,
                             CandidateRegions* parent_cr) 
{
    NECNode* leaf = path->back();
    ULONG child_score = MatchingOrderCRScore(child_cr, nontree_edges, leaf);
    ULONG parent_score = MatchingOrderCRScore(parent_cr, nontree_edges, leaf);
#ifdef PRINT_ON
    cout << "Matching Order Score " << child_score + parent_score << endl;
#endif
    return child_score + parent_score;
}

vector<NECNode*>* GenerateNECOrder(MatchingOrderPq* path_pq, 
                                   vector<NECNode*>& initialized_nec) 
{
    unordered_set<NECNode*> seen_nec;
    vector<NECNode*>* nec_queue;
    MatchingOrderPair* pair;
    NECNode* curr_nec;
    vector<NECNode*>* evaluation_queue = new vector<NECNode*>;

    for (int i = 0; i < initialized_nec.size(); i += 1) {
        seen_nec.insert(initialized_nec[i]);
        evaluation_queue->push_back(initialized_nec[i]);
    }

    while (!path_pq->empty()) {
        pair = path_pq->top();
        path_pq->pop();
        nec_queue = pair->q;
        
        for (int index = nec_queue->size() - 1; index >= 0; index -= 1) {
            curr_nec = nec_queue->at(index);
            if (seen_nec.count(curr_nec) > 0 || curr_nec->members->size() == 0) {
                continue;
            }

            seen_nec.insert(curr_nec);
            evaluation_queue->push_back(curr_nec);
        }
    }

    return evaluation_queue;
}

void MatchingOrderPrint(MatchingOrderPq* pq) {
    cout << " Matching Order " << endl;
    cout << "Pq Size " << pq->size() << endl;
    while (!pq->empty()) {
        MatchingOrderPair* pair = pq->top();
        pq->pop();
        vector<NECNode*>* q = pair->q;

        cout << "path size " << q->size() << endl;
        for (int i = 0; i < q->size(); i += 1) {
            cout << pair << " " << q->at(i) << " - ";
        }
        cout << endl;
    }
}
bool MatchingOrderNECCompare(pair<NECNode*, int>& p1, pair<NECNode*, int>& p2) {
    return p1.second < p2.second;
}

vector<NECNode*>* MatchingOrderByNEC(NECNode* nec_root, CandidateRegions* cr) {
    NECNode* nec;
    vector<NECNode*>* nec_queue = new vector<NECNode*>;
    vector<pair<NECNode*, int> > sort_buffer;
    unordered_map<VERTEX,vector<VERTEX>*>* regions;
    unordered_map<VERTEX,vector<VERTEX>*>::iterator region_it;
    int score;
    CandidateRegions::iterator it = cr->begin();
    while (it != cr->end()) {
        nec = it->first;
        regions = it->second;
        region_it = regions->begin();
        score = 0;
        while (region_it != regions->end()) {
            score += region_it->second->size();
            region_it++;
        }
        sort_buffer.push_back(make_pair<NECNode*, int>(nec, score));
        ++it;
    }

    sort(sort_buffer.begin(), sort_buffer.end(), MatchingOrderNECCompare);
    for (int i = 0; i < sort_buffer.size(); i += 1) {
        nec_queue->push_back(sort_buffer[i].first);
    }
    return nec_queue;
}



vector<NECNode*>* MatchingOrderByPath(NECNode* nec_root, 
                                      CandidateRegions* child_cr,
                                      CandidateRegions* parent_cr,
                                      vector<NECNode*>& initialized_nec) 
{
    ULONG new_nontree_edges;
    ULONG nontree_edges;
    vector<NECNode*>* path;
    NECNode* leaf_node;
    NECNode* child;
    vector<NECNode*>* NEC_to_visit;
    pair<vector<NECNode*>*, ULONG>* p;
    pair<vector<NECNode*>*, ULONG>* new_pair;
    MatchingOrderPair* pq_elem;

    queue<pair<vector<NECNode*>*, ULONG>* > bfs_queue;
    MatchingOrderPq* pq = new MatchingOrderPq;
    unordered_set<NECNode*> visited;
    vector<NECNode*>* seed_queue = new vector<NECNode*>;

    new_pair = new pair<vector<NECNode*>*, ULONG>;
    new_pair->first = seed_queue;
    new_pair->second = 0;
    seed_queue->push_back(nec_root);
    bfs_queue.push(new_pair);
    visited.insert(nec_root);

    while (bfs_queue.empty() == false) {
        p = bfs_queue.front();
        bfs_queue.pop();
        path = p->first;
        leaf_node = path->back();
        NEC_to_visit = MatchingOrderGetUnvisited(leaf_node, &visited);
        new_nontree_edges = leaf_node->parents->size() +
                            leaf_node->children->size() -
                            NEC_to_visit->size();
        nontree_edges = p->second + new_nontree_edges;

        for (int i = 1; i < NEC_to_visit->size(); i += 1) {
            child = NEC_to_visit->at(i);
            MatchingOrderVisitChild(child, path, nontree_edges, &bfs_queue, &visited);
        }

        if (NEC_to_visit->size() > 0) {
            // Reuse allocated queue for 0-th child.
            child = NEC_to_visit->at(0);
            MatchingOrderVisitChild(child, path, nontree_edges, &bfs_queue, &visited);
        }

        if (NEC_to_visit->size() == 0) {
            // End of path
            pq_elem = new MatchingOrderPair;
            pq_elem->q = path;
            pq_elem->score = MatchingOrderCalcScore(path, nontree_edges, 
                                                    child_cr, parent_cr);
            pq->push(pq_elem);
        }
    }

    return GenerateNECOrder(pq, initialized_nec);
}

void PrintCandidateRegions(CandidateRegions* cr) {
    CandidateRegions::iterator it = cr->begin();
    while (it != cr->end()) {
        cout << "NEC NODE :" << it->first << " - \n";
        unordered_map<VERTEX, vector<VERTEX>*>::iterator it2 = it->second->begin();
        while (it2 != it->second->end()) {
            cout << "\t\tDataNode vertex prior: " << it2->first << " - ";
            vector<VERTEX>* v = it2->second;
            for (int i = 0; i < v->size(); i += 1) {
                cout << v->at(i) << " ";
            }
            cout << endl;
            it2++;
        }
        it ++;
    }
}

/*
 * Parameters:
 * dg, qg :
 *      the data graph and query graph data types, respectively.
 *
 * nec_root :
 *      the starting NEC node to generate candidate regions
 *
 * child_cr and parent_cr :
 *      data structures for storing the candidate
 *      regions for initialized parent and children nec vertices, respectively.
 *
 * initialized_nec : 
 *      a collection of nec_vertices that have already been
 *      initialized with a data vertex
 */
void AllocCandidateRegions(DataGraph& dg,
                            QueryGraph* qg,
                            NECNode* nec_root,
                            CandidateRegions*& child_cr,
                            CandidateRegions*& parent_cr,
                            vector<NECNode*>& initialized_nec) 
{
    unordered_set<VERTEX>* visited_vertices = new unordered_set<VERTEX>;
    unordered_set<NECNode*> visited_nec;
    LABEL root_label;
    VERTEX root_node_id;
    vector<VERTEX> candidate_data_vertices;
    QueryNode* allocated_query_node;
    VERTEX mapped_node;
    NECNode* child_nec;
    VERTEX* data_node_iter;
    bool initialized_query_node_has_cr = true;

    if (nec_root->members->size() == 0) {
        /*
         * Handle the initialized query vertex case where the NEC root is an
         * empty supernode with children NEC containing singleton query vertices
         * that have been initialized
         */
        for (int child = 0; child < nec_root->children->size(); child += 1) {
            child_nec = nec_root->children->at(child);
            allocated_query_node = child_nec->members->at(0);
            candidate_data_vertices.clear();
            mapped_node = qg->initialized_query_nodes->at(allocated_query_node);
            candidate_data_vertices.push_back(mapped_node);

            // NEC Nodes do keep parent edges to the supernode. This is a
            // special case so graph traversal does not ever visit the
            // meaningless supernode which does not correspond to any query
            // vertex.
            if (!ExploreCandidateRegions(dg, qg, child_nec, 
                                         &candidate_data_vertices,
                                         visited_vertices, child_cr, parent_cr, 
                                         child_cr, -1, NULL, &visited_nec))
            {
                initialized_query_node_has_cr = false;
                break;
            }
        }
        if (!initialized_query_node_has_cr) {
            // DeleteCandidateRegion(child_cr);
            // DeleteCandidateRegion(parent_cr);
            child_cr = NULL;
            parent_cr = NULL;
        }
    } else {
        root_label = nec_root->members->at(0)->label;
        data_node_iter = dg.vfs->GetVertexIterator(root_label);
        for(int data_node_index = 0; 
            data_node_index < dg.vfs->NumVertices(root_label); 
            data_node_index += 1)
        {
            // TODO : FILTER data nodes based on NEC characteristics.
            root_node_id = *(data_node_iter + data_node_index);
            candidate_data_vertices.clear();
            candidate_data_vertices.push_back(root_node_id);

            if(!ExploreCandidateRegions(dg, qg, nec_root,
                                        &candidate_data_vertices,
                                        visited_vertices, child_cr,
                                        parent_cr, child_cr,
                                        -1, NULL, &visited_nec))
            {
                continue;
            }
        }
    }
}

// Verify that a candidate data node has all outgoing edges to the already
// matched parnt nodes.
bool SubGraphParentReferences(vector<VERTEX>* parents, 
                                VERTEX data_node,
                                DataGraph& dg) {
    VERTEX parent;
    for (int i = 0; i < parents->size(); i += 1) {
        parent = parents->at(i);
        if (!dg.efs->HasEdge(parent, data_node)) {
            return false;
        }
    }
    return true;
}

// Verify that a candidate data node has all outgoing edges to the already
// matched children nodes.
bool SubGraphChildrenReferences(vector<VERTEX>* children, 
                                VERTEX data_node,
                                DataGraph& dg) {
    VERTEX child;
    for (int i = 0; i < children->size(); i += 1) {
        child = children->at(i);
        if (!dg.efs->HasEdge(data_node, child)) {
            return false;
        }
    }
    return true;
}

vector<VERTEX>* SubGraphDataVertexUnion(vector<vector<VERTEX>*>* collections,
                                        unordered_set<VERTEX>* marked_data_vertex,
                                        vector<VERTEX>* matched_child_dnodes,
                                        vector<VERTEX>* matched_parent_dnodes,
                                        DataGraph& dg) 
{
    vector<VERTEX>* vec;
    VERTEX data_node;
    VERTEX matched_node;
    unordered_set<VERTEX>* result = new unordered_set<VERTEX>;
    for (int i = 0; i < collections->size(); i += 1) {
        vec = collections->at(i);
        for (int j = 0; j < vec->size(); j += 1) {
            data_node = vec->at(j);
            if (marked_data_vertex->count(data_node) == 0 &&
                SubGraphChildrenReferences(matched_child_dnodes, data_node, dg) &&
                SubGraphParentReferences(matched_parent_dnodes, data_node, dg))
            {
                result->insert(data_node);
            }
        }
    }

    for (int i = 0; i < matched_child_dnodes->size(); i += 1) {
        matched_node = matched_child_dnodes->at(i);
        if (result->count(matched_node) > 0) {
            result->erase(matched_node);
        }
    }
    vector<VERTEX>* ret = new vector<VERTEX>(result->begin(), result->end());
    return ret;
}

unordered_set<VERTEX>* SubGraphDataVertexIntersection(vector<vector<VERTEX>*>* collections,
                                               unordered_set<VERTEX>* marked_data_vertex,
                                               vector<VERTEX>* matched_child_dnodes,
                                               vector<VERTEX>* matched_parent_dnodes,
                                               DataGraph& dg) 
{
    if (collections->size() == 0)
        return new unordered_set<VERTEX>;

    vector<VERTEX>* seed = collections->at(0);
    VERTEX matched_node;
    vector<VERTEX>* ret;
    unordered_set<VERTEX>* result = new unordered_set<VERTEX>();
    unordered_set<VERTEX>::iterator it;

    for (int i = 0; i < seed->size(); i += 1) {
        matched_node = seed->at(i);
        if (marked_data_vertex->count(matched_node) == 0 &&
            SubGraphChildrenReferences(matched_child_dnodes, matched_node, dg) &&
            SubGraphParentReferences(matched_parent_dnodes, matched_node, dg))
        {
            result->insert(seed->at(i));
        }
    }

    for (int i = 1; i < collections->size(); i += 1) {
        seed = collections->at(i);
        unordered_set<VERTEX> cmp = unordered_set<VERTEX>(seed->begin(), seed->end());
        it = result->begin();
        while (it != result->end()) {
            if (cmp.count(*(it)) == 0) {
                it = result->erase(it);
            } else {
                it++;
            }
        }
    }
    return result;
}

// Sort datanode vertices for std::set_intersection
void SortCr(CandidateRegions* cr) {
    CandidateRegions::iterator it1 = cr->begin();
    unordered_map<VERTEX,vector<VERTEX>*>::iterator it2;

    while (it1 != cr->end()) {
        it2 = it1->second->begin();
        while (it2 != it1->second->end()) {
            sort(it2->second->begin(), it2->second->end());
        }
    }
}

/*
 * Returns a vector of all data vertices mapped to query vertices of
 * all neighbor necs
 */
vector<VERTEX>* SubGraphGetMatchedChildren(NECNode* nec, 
                                           unordered_set<NECNode*>* matched_nec,
                                           unordered_map<QueryNode*, VERTEX>* mapping) 
{
    vector<VERTEX>* matched_child_dnodes = new vector<VERTEX>;
    NECNode* neighbor;
    QueryNode* query_node;

    // Get matched child query vertices
    for (int i = 0; i < nec->children->size(); i += 1) {
        neighbor = nec->children->at(i);
        if (matched_nec->count(neighbor) == 0)
            continue;

        for (int member_index = 0;
             member_index < neighbor->members->size();
             member_index += 1)
        {
            query_node = neighbor->members->at(member_index);
            matched_child_dnodes->push_back(mapping->at(query_node));
        }
    }
    return matched_child_dnodes;
}

vector<VERTEX>* SubGraphGetMatchedParents(NECNode* nec, 
                                          unordered_set<NECNode*>* matched_nec,
                                          unordered_map<QueryNode*, VERTEX>* mapping) 
{
    vector<VERTEX>* matched_parent_dnodes = new vector<VERTEX>;
    NECNode* neighbor;
    QueryNode* query_node;

    // Get matched child query vertices
    for (int i = 0; i < nec->parents->size(); i += 1) {
        neighbor = nec->parents->at(i);
        if (matched_nec->count(neighbor) == 0)
            continue;

        for (int member_index = 0;
             member_index < neighbor->members->size();
             member_index += 1)
        {
            query_node = neighbor->members->at(member_index);
            matched_parent_dnodes->push_back(mapping->at(query_node));
        }
    }
    return matched_parent_dnodes;
}

vector<vector<VERTEX>*>* CollectCandidateRegions(NECNode* nec, 
                                                 CandidateRegions* cr)
{
    SubRegions* NEC_regions;
    SubRegions::iterator subregion_it;
    vector<vector<VERTEX>*>* result = new vector<vector<VERTEX>*>;
    if (cr->count(nec) == 0)
        return result;
    NEC_regions = cr->at(nec);
    subregion_it = NEC_regions->begin();
    while (subregion_it != NEC_regions->end()) {
        result->push_back(subregion_it->second);
        subregion_it++;
    }
    return result;
}

bool HasEdgesTo(VERTEX v, 
                DataGraph& dg, 
                vector<VERTEX>* matched_children_dnodes,
                vector<VERTEX>* matched_parent_dnodes)
{
    for (int i = 0; i < matched_children_dnodes->size(); i += 1) {
        if (!dg.efs->HasEdge(v, matched_children_dnodes->at(i)))
            return false;
    }
    for (int i = 0; i < matched_parent_dnodes->size(); i += 1) {
        if (!dg.efs->HasEdge(matched_parent_dnodes->at(i), v))
            return false;
    }
    return true;
}

unordered_set<VERTEX>* FlattenCandidateRegions(NECNode* nec,
                                               CandidateRegions* cr,
                                               vector<VERTEX>* matched_children_dnodes,
                                               vector<VERTEX>* matched_parent_dnodes,
                                               DataGraph& dg)
{
    unordered_set<VERTEX>* result = new unordered_set<VERTEX>;
    SubRegions* NEC_regions;
    vector<VERTEX>* candidates;
    SubRegions::iterator subregion_it;
    if (cr->count(nec) == 0)
        return result;
    NEC_regions = cr->at(nec);
    subregion_it = NEC_regions->begin();
    while (subregion_it != NEC_regions->end()) {
        candidates = subregion_it->second;
        for (int i = 0; i < candidates->size(); i += 1) {
            if (HasEdgesTo(candidates->at(i), dg,
                           matched_children_dnodes, 
                           matched_parent_dnodes))
                result->insert(candidates->at(i));
        }
        subregion_it++;
    }
    return result;
}

bool SubGraphCountCRCmp(vector<VERTEX>* v1, vector<VERTEX>* v2) {
    return (v1->size() < v2->size());
}

unordered_set<VERTEX>* SetIntersection(unordered_set<VERTEX>* s1, 
                                       unordered_set<VERTEX>* s2)
{
    unordered_set<VERTEX>::iterator it = s1->begin();
    while (it != s1->end()) {
        if (s2->count( *(it) ) == 0) {
            it = s1->erase(it);
        } else {
            it ++;
        }
    }
    return s1;
}

/*
 * 1. Intersect all child CR for the nec
 * 2. Intersect all parent CR for the nec
 * 3. Intersect 1 and 2.
 */
vector<VERTEX>* SubGraphSmallestCRIntersection(NECNode* nec,
                                    unordered_set<VERTEX>* marked_data_vertex,
                                    CandidateRegions* child_cr,
                                    CandidateRegions* parent_cr,
                                    DataGraph& dg)
{
    vector<VERTEX> matched_children;
    vector<VERTEX> matched_parents;
    unordered_set<VERTEX>* data_nodes;
    vector<VERTEX>* result;
    vector<vector<VERTEX>*>* child_candidates = CollectCandidateRegions(nec, child_cr);
    vector<vector<VERTEX>*>* parent_candidates = CollectCandidateRegions(nec, parent_cr);
    sort(child_candidates->begin(), child_candidates->end(), SubGraphCountCRCmp);
    sort(parent_candidates->begin(), parent_candidates->end(), SubGraphCountCRCmp);
    unordered_set<VERTEX>* child_intersection = 
                            SubGraphDataVertexIntersection(child_candidates,
                                                           marked_data_vertex,
                                                           &matched_children,
                                                           &matched_parents,
                                                           dg);
    unordered_set<VERTEX>* parent_intersection = 
                            SubGraphDataVertexIntersection(parent_candidates,
                                                           marked_data_vertex,
                                                           &matched_children,
                                                           &matched_parents,
                                                           dg);
 
    if (parent_intersection->size() < child_intersection->size())
        data_nodes = SetIntersection(parent_intersection, child_intersection);
    else
        data_nodes = SetIntersection(child_intersection, parent_intersection);

    result = new vector<VERTEX>(data_nodes->begin(), data_nodes->end());
    delete child_candidates;
    delete parent_candidates;
    delete child_intersection;
    delete parent_intersection;
    return result;
}

/*
 * Gets candidate nodes for NEC if neighbor NEC have been matched already by
 * enforcing all candidate datanodes to have edges to the matched data vertices.
 */
vector<VERTEX>* SubGraphCandidatesWithPriors(NECNode* nec,
                                    unordered_set<VERTEX>* marked_data_vertex,
                                    CandidateRegions* child_cr,
                                    CandidateRegions* parent_cr,
                                    vector<VERTEX>* matched_children_dnodes,
                                    vector<VERTEX>* matched_parent_dnodes,
                                    DataGraph& dg)
{
    unordered_set<VERTEX>* child_candidates = 
                        FlattenCandidateRegions(nec, 
                                                child_cr, 
                                                matched_children_dnodes,
                                                matched_parent_dnodes,
                                                dg);
    unordered_set<VERTEX>* parent_candidates = 
                        FlattenCandidateRegions(nec,
                                                parent_cr,
                                                matched_children_dnodes,
                                                matched_parent_dnodes, 
                                                dg);
    unordered_set<VERTEX>* intersect;
    vector<VERTEX>* result;
    if (child_candidates->size() < parent_candidates->size())
        intersect = SetIntersection(child_candidates, parent_candidates);
    else
        intersect = SetIntersection(parent_candidates, child_candidates);
    result = new vector<VERTEX>(intersect->begin(), intersect->end());
    delete intersect;
    delete child_candidates;
    delete parent_candidates;
    return result;
}

// Get candidate data vertex for a NEC node
vector<VERTEX>* SubGraphCountGetDataVertex(NECNode* nec, 
                                unordered_map<QueryNode*, VERTEX>* mapping,
                                unordered_set<NECNode*>* matched_nec,
                                unordered_set<VERTEX>* marked_data_vertex,
                                CandidateRegions* child_cr,
                                CandidateRegions* parent_cr,
                                DataGraph& dg)
{
    NECNode* neighbor;
    QueryNode* query_node;
    VERTEX mapped_data_node;
    unordered_map<VERTEX, vector<VERTEX>*>* nec_cr_map;
    unordered_map<VERTEX, vector<VERTEX>*>::iterator cr_iter;
    vector<vector<VERTEX>*> valid_candidate_regions;
    vector<VERTEX>* matched_parent_dnodes;
    vector<VERTEX>* matched_child_dnodes;
    vector<VERTEX>* result;

    // NEC doesn't exist. Means there was no valid CR.
    /*
    if (cr->count(nec) == 0) {
        return new vector<VERTEX>;
    }
    */

    matched_child_dnodes = SubGraphGetMatchedChildren(nec, matched_nec, mapping);
    matched_parent_dnodes = SubGraphGetMatchedParents(nec, matched_nec, mapping);
    int num_matched = 0;

    if (matched_child_dnodes->size() == 0 && matched_parent_dnodes->size() == 0)
    {
        result = SubGraphSmallestCRIntersection(nec, marked_data_vertex, child_cr,
                                                parent_cr, dg);
    } else {
        result = SubGraphCandidatesWithPriors(nec, marked_data_vertex, child_cr,
                                              parent_cr, matched_child_dnodes,
                                              matched_parent_dnodes, dg);
    }
    return result;
}

bool SubGraphIsJoinable(QueryNode* qnode, VERTEX dnode, 
                        unordered_map<QueryNode*, VERTEX>* mapping,
                        DataGraph& dg) 
{
    QueryNode* adj_node;
    VERTEX mapped_node;

    // Verify parents. Most likely neglected edges from NEC construction.
    for (int p_index = 0; p_index < qnode->parents->size(); p_index += 1) {
        adj_node = qnode->parents->at(p_index);
        if (mapping->count(adj_node) > 0) {
            mapped_node = mapping->at(adj_node);
            if (!dg.efs->HasEdge(mapped_node, dnode)) {
                return false;
            }
        }
    }

    for (int c_index = 0; c_index < qnode->children->size(); c_index += 1) {
        adj_node = qnode->children->at(c_index);
        if (mapping->count(adj_node) > 0) {
            mapped_node = mapping->at(adj_node);
            if (!dg.efs->HasEdge(dnode, mapped_node)) {
                return false;
            }
        }
    }
    return true;
}

void SubGraphMapping(unordered_map<QueryNode*, VERTEX>* mapping,
                     vector<VERTEX>* combination,
                     NECNode* nec_node,
                     unordered_set<NECNode*>* matched_nec,
                     unordered_set<VERTEX>* marked_data_vertex)
{
    QueryNode* qnode;
    matched_nec->insert(nec_node);
    for (int i = 0; i < nec_node->members->size(); i += 1) {
        qnode = nec_node->members->at(i);
        mapping->insert(make_pair<QueryNode*, VERTEX>(qnode, combination->at(i)));
        marked_data_vertex->insert(combination->at(i));
    }

}
void SubGraphRestore(unordered_map<QueryNode*, VERTEX>* mapping,
                     vector<VERTEX>* combination,
                     NECNode* nec_node,
                     unordered_set<NECNode*>* matched_nec,
                     unordered_set<VERTEX>* marked_data_vertex)
{
    QueryNode* qnode;
    matched_nec->erase(nec_node);
    for (int i = 0; i < nec_node->members->size(); i += 1) {
        qnode = nec_node->members->at(i);
        mapping->erase(qnode);
        marked_data_vertex->erase(combination->at(i));
    }

}

void SubGraphPrintMapping(unordered_map<QueryNode*, VERTEX>* mapping) {
    unordered_map<QueryNode*, VERTEX>::iterator it = mapping->begin();
    while (it != mapping->end()) {
        cout << it->second << " - ";
        it++;
    }
    cout << endl;
}

void SubGraphPrintMapping(const string output_file,
                          unordered_map<QueryNode*, VERTEX>* mapping) 
{
    ofstream outfile;
    outfile.open(output_file, ios_base::app);

    unordered_map<QueryNode*, VERTEX>::iterator it = mapping->begin();
    while (it != mapping->end()) {
        outfile << it->first << " - " << it->second << endl;
        it++;
    }
    outfile << endl;
}

int SubGraphCount(vector<NECNode*>* matching_order, int matching_index,
                   unordered_map<QueryNode*, VERTEX>* mapping, 
                   unordered_set<NECNode*>* matched_nec,
                   unordered_set<VERTEX>* marked_data_vertex,
                   CandidateRegions* child_cr,
                   CandidateRegions* parent_cr,
                   DataGraph& dg) 
{
    NECNode* curr_nec;
    vector<VERTEX>* candidate_data_vertex;
    vector<VERTEX>* combination;
    bool matched;
    int count = 0;
    curr_nec = matching_order->at(matching_index);

    // Intersection of (parent_cr <NEC,priors> union child_cr <NEC,priors>)
    candidate_data_vertex = SubGraphCountGetDataVertex(curr_nec, mapping,
                                                       matched_nec,
                                                       marked_data_vertex,
                                                       child_cr,
                                                       parent_cr,
                                                       dg);

#ifdef PRINT_ON
    if (candidate_data_vertex->size() == 0) {
        cout << "no candidate data vertex for nec :" << curr_nec << endl;
    }
    cout << endl;
    cout << "NEC:" << curr_nec << " CANIDATE SIZE" << candidate_data_vertex->size() <<  endl;
    cout << "DATA VERTEX CANDIDATE:";
    for (int i = 0; i < candidate_data_vertex->size(); i += 1) {
        cout << candidate_data_vertex->at(i) << " - ";
    }
    cout << endl;
#endif
    XChoseY combinations = XChoseY(candidate_data_vertex, curr_nec->members->size());
    combination = combinations.getNext();

    while (combination != NULL) {
        matched = true;
        for (int i = 0; i < combination->size(); i += 1) 
        {
            if (!SubGraphIsJoinable(curr_nec->members->at(i), 
                                   combination->at(i),
                                   mapping,
                                   dg))
            {
                matched = false;
                break;
            }
        }

        if (!matched) {
            combination = combinations.getNext();
            continue;
        }

        SubGraphMapping(mapping, combination, curr_nec, matched_nec, 
                        marked_data_vertex);

        if (matching_index == matching_order->size() - 1) {
            //SubGraphPrintMapping(mapping);
#if defined PRINT_SUBGRAPH_MATCHED && defined SUBGRAPH_OUTPUT_FILE
            SubGraphPrintMapping(SUBGRAPH_OUTPUT_FILE, mapping);
#endif
            count += 1;
        } else {
            count += SubGraphCount(matching_order, matching_index + 1,
                                   mapping, matched_nec,
                                   marked_data_vertex, child_cr, parent_cr,
                                   dg);
        }
        SubGraphRestore(mapping, combination, curr_nec, matched_nec,
                        marked_data_vertex);

        combination = combinations.getNext();
    }
    return count;
}

void PrintNECTreeH(NECNode* root, unordered_set<NECNode*>& seen) {
    NECNode* curr;
    seen.insert(root);
    if (root == NULL)
        return;
    cout << "NEC NODE: " << root;
    cout << " Members: ";
    for (int i = 0; i < root->members->size(); i += 1) {
        cout << root->members->at(i) << "-" ;
    }
    cout << " Parent NEC: ";
    for (int i = 0; i < root->parents->size(); i += 1) {
        cout << root->parents->at(i) << "-";
    }
    cout << " Children NEC: ";
    for (int i = 0; i < root->children->size(); i += 1) {
        cout << root->children->at(i) << "-";
    }
    cout << endl;

    for (int i = 0; i < root->children->size(); i += 1) {
        curr = root->children->at(i);
        if (seen.count(curr) > 0)
            continue;
        seen.insert(curr);
        PrintNECTreeH(curr, seen);
    }
    for (int i = 0; i < root->parents->size(); i += 1) {
        curr = root->parents->at(i);
        if (seen.count(curr) > 0)
            continue;
        seen.insert(curr);
        PrintNECTreeH(curr, seen);
    }
}

void PrintNECTree(NECNode* root) {
    unordered_set<NECNode*> seen;
    PrintNECTreeH(root, seen);
}

int TurboIso(DataGraph& dg, QueryGraph* query_graph) {
    CandidateRegions* child_cr = new CandidateRegions;
    CandidateRegions* parent_cr = new CandidateRegions;
    query_node_map* qg = query_graph->query_nodes;
    QueryNode* start_vertex;
    LABEL root_label;
    NECNode* nec_root;
    MatchingOrderPq* matching_order_pq;
    vector<NECNode*>* matching_queue;
    unordered_set<NECNode*> matched_NEC;
    unordered_map<QueryNode*, VERTEX> query_node_mapping;
    unordered_set<VERTEX> marked_data_vertex;
    vector<NECNode*> initialized_nec;
    int result = 0;

    start_vertex = GetStartingQueryVertex(dg, query_graph);
    nec_root = RewriteToNECTree(query_graph, start_vertex, initialized_nec);

    PrintNECTree(nec_root);

    // Takes argument of matched query nodes-NEC and generates special CR.
    AllocCandidateRegions(dg, query_graph, nec_root, child_cr,
                          parent_cr, initialized_nec);
    if (child_cr != NULL && parent_cr != NULL) {
#ifdef PRINT_ON
        cout << " Child CR " << endl;
        PrintCandidateRegions(child_cr);
        cout << " Parent CR " << endl;
        PrintCandidateRegions(parent_cr);
#endif
        // Matching queue needs to look at both candidate regions.
        matching_queue = MatchingOrderByPath(nec_root, child_cr,
                                             parent_cr, initialized_nec);
        //matching_queue = MatchingOrderByNEC(nec_root, cr);

        matched_NEC.clear();
        query_node_mapping.clear();
        // SubGraphCount cannot treat child and parent vertices using the same
        // candidate region
        result = SubGraphCount(matching_queue, 0, &query_node_mapping, 
                               &matched_NEC, &marked_data_vertex, 
                               child_cr, parent_cr, dg);
    }

    if (child_cr == NULL || parent_cr == NULL)
        cout << "cr null" << endl;

    return result;
}

void QueryAllAuthors(string author_list, string query_graph_file, string output_file, DataGraph& dg) {
    int result;
    QueryGraph* qg5;
    string line;
    clock_t clock_begin;
    clock_t clock_end;
    double elapsed_time;
    qg5 = ReadQueryGraphFromFile(query_graph_file);
    ifstream author_file(author_list);
    ofstream querygraph5(output_file);

    while (getline(author_file, line)) {
        unordered_map<QueryNode*, UINT>::iterator it;
        it = qg5->initialized_query_nodes->begin();
        it->second = stoul(line);

        clock_begin = clock();
        result = TurboIso(dg, qg5);
        clock_end = clock();
        elapsed_time = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
        querygraph5 << 1000*elapsed_time << "," << result << endl;

        querygraph5.flush();
    }
    querygraph5.close();
}

/*
 * Runs algorithm on data graph and query graph specified, including the
 * initialized query nodes in the query_graph_file specification
 */
void QuerySingleInstance(string query_graph_file, string output_file, DataGraph& dg) {
    int result;
    QueryGraph* qg5;
    clock_t clock_begin;
    clock_t clock_end;
    double elapsed_time;
    qg5 = ReadQueryGraphFromFile(query_graph_file);
    ofstream querygraph5(output_file);

    //cout << qg5->query_nodes->size() << "size" << endl;
    PrintQueryGraph(qg5);

    clock_begin = clock();
    result = TurboIso(dg, qg5);
    clock_end = clock();
    elapsed_time = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
    querygraph5 << 1000*elapsed_time << "," << result << endl;
    querygraph5.close();
}

int main(int argc, char** argv) {
    DataGraph dg;
    QueryGraph* qg5;
    string line;
    int result;
    clock_t clock_begin;
    clock_t clock_end;
    double elapsed_time;

    /*
    if (argc < 1)
        // Need at least query graph
        return 1;
    */

    /*
    VertexFs vfs("dblp_vfs");
    EdgeFs efs("dblp_efs.txt", "dblp_efs_rev.txt", &vfs);
    */
    VertexFs vfs("data_graphs/vertex_file");
    EdgeFs efs("data_graphs/edge_file", "data_graphs/reverse_edge_file", &vfs);
    
    dg.vfs = &vfs;
    dg.efs = &efs;

    //qg = ReadQueryGraphFromFile(argv[1]);
    QuerySingleInstance("query_graphs/yloop.txt", "output/yloop_single", dg);
    //QueryAllAuthors("author_list", "query_graphs/query_graph2.txt", "output/qg2_out", dg);
}
#endif
