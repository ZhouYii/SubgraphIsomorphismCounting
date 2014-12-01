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
bool ExploreCandidateRegions(DataGraph& dg, NECNode* nec_root,
                             vector<VERTEX>* data_node_candidates,
                             unordered_set<VERTEX>* marked_vertex, 
                             CandidateRegions* cr,
                             VERTEX parent);

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

/*
void GetSortedNeighbors(vector<QueryNode*>* nodes, vector<QueryNode*>* els) {
    for(int i = 0; i < nodes->size(); i += 1)
        labels->push_back(nodes->at(i));

    // Sort Ascending. Order doesn't actually matter.
    sort(labels->begin(), labels->end());
}
*/

vector<NECNode*>* MatchNECNodes(NECBucketMap* node_buckets) {
    vector<NECNode*>* nec_nodes;
    vector<QueryNode*>* query_nodes;
    vector<LABEL> nec_child_labels;
    vector<LABEL> nec_parent_labels;
    vector<LABEL> cmp_child_labels;
    vector<LABEL> cmp_parent_labels;
    QueryNode* query_node;
    QueryNode* origin_node;
    QueryNode* cmp_node;
    NECNode* nec_node;

    nec_nodes = new vector<NECNode*>;
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
            nec_nodes->push_back(nec_node);

            // Get adjacent information
            // This is the slow path. Take this after the hash path :
            // TODO : two separate hashes for parent and child labels

            // Iterate over all remaining nodes in vector (in reverse) to find
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

    return nec_nodes;
}

vector<NECNode*>* GenerateNECNodes(unordered_set<QueryNode*>* pnodes,
                                   unordered_set<QueryNode*>* seen)
{
    vector<QueryNode*>* bucket_vector;
    NECBucketMap query_node_buckets;
    QueryNode* curr_node;
    
    unordered_set<QueryNode*>::iterator it = pnodes->begin();
    // Group Query Nodes by NEC
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

NECNode* RewriteToNECTree(query_node_map* qg, QueryNode* start) {
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

            /*
            cout << "NEC Node ID :" << curr << " Children:";
            for (int i = 0; i < curr->children->size(); i += 1) {
                cout << curr->children->at(i) << " | ";
            }
            cout << " \t\tParents:";
            for (int i = 0; i < curr->parents->size(); i += 1) {
                cout << curr->parents->at(i) << " | ";
            }
            cout << " \t\t Members:";
            for (int i = 0; i < curr->members->size(); i += 1) {
                cout << curr->members->at(i) << " | ";
            }
            cout << endl;
            */
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

NEC_PQ* GetNECEvalOrder(DataGraph& dg, NECNode* nec_root, 
                        VERTEX v_index, int parent_flag) {
    vector<NECNode*>* nec_nodes;
    vector<NECNode*>::iterator it;
    pair<NECNode*, UINT>* p;
    NECNode* curr;
    NEC_PQ* pq = new NEC_PQ;

    if (parent_flag)
        nec_nodes = nec_root->parents;
    else 
        nec_nodes = nec_root->children;
    for (it = nec_nodes->begin(); it != nec_nodes->end(); it++) {
        curr = *(it);
        p = new pair<NECNode*, UINT>;
        p->first = curr;
        if (parent_flag)
            p->second = dg.efs->GetInDegree(v_index, curr->members->at(0)->label);
        else
            p->second = dg.efs->GetOutDegree(v_index, curr->members->at(0)->label);
        pq->push(p);
    }

    return pq;
}

void ClearCandidateRegions(CandidateRegions* cr, 
                           vector<NECNode*>* visited_nec_neighbors,
                           VERTEX parent)
{
    NECNode* clear_nec;

    for (int neighbor = 0; neighbor < visited_nec_neighbors->size(); neighbor += 1) 
    {
        clear_nec = visited_nec_neighbors->at(neighbor);
        if (cr->count(clear_nec) > 0 &&
            cr->find(clear_nec)->second->count(parent) > 0) {
            cr->find(clear_nec)->second->find(parent)->second->clear();
        }
    }
}

void InsertIntoCandidateRegion(CandidateRegions* cr, 
                               VERTEX parent,
                               NECNode* nec_node,
                               VERTEX vertex_id) {
    unordered_map<VERTEX, vector<VERTEX>*>* prior_vertex_map;

    //cout << "Insert into candidate regions: parent - " << parent << " vertex - " << nec_node << endl;

    if (cr->count(nec_node) == 0) {
        cr->insert(make_pair<NECNode*, unordered_map<VERTEX, vector<VERTEX>*>*>
                        (nec_node, new unordered_map<VERTEX, vector<VERTEX>*>));
        //cout << "CR ADD " << nec_node << endl;
    }

    prior_vertex_map = cr->find(nec_node)->second;
    if (prior_vertex_map->count(parent) == 0) {
        prior_vertex_map->insert(make_pair<VERTEX, vector<VERTEX>*>
                                          (parent, new vector<VERTEX>));
        //cout << "CR PARENT ADD " << parent << endl;
    }
    prior_vertex_map->find(parent)->second->push_back(vertex_id);
}

bool ExploreSubRegions(DataGraph& dg,
                       unordered_set<VERTEX>* marked_vertex, 
                       CandidateRegions* cr,
                       VERTEX parent,
                       VERTEX adj_vertex,
                       vector<NECNode*>* visited_nec_neighbors,
                       NEC_PQ* pq)
{
    NECNode* adj_nec;
    vector<UINT>* candidate_vertices;
    LABEL candidate_label;
    bool matched = true;

    if (pq->size() > 0) {
        adj_nec = pq->top()->first;
        candidate_label = adj_nec->members->at(0)->label;
        candidate_vertices = dg.efs->GetOutdegreeVertices(adj_vertex, 
                                                        candidate_label);
        // Construct candidate vertex using VFS
        if (ExploreCandidateRegions(dg, adj_nec, candidate_vertices,
                                    marked_vertex,
                                    cr, adj_vertex) == false)
        {
            ClearCandidateRegions(cr, visited_nec_neighbors, parent);
            matched = false;
        }
        visited_nec_neighbors->push_back(adj_nec);
        delete pq->top();
        pq->pop();
        delete candidate_vertices;
    }
    return matched;
}

bool ExploreCandidateRegions(DataGraph& dg, NECNode* nec_root,
                             vector<VERTEX>* data_node_candidates,
                             unordered_set<VERTEX>* marked_vertex, 
                             CandidateRegions* cr,
                             VERTEX parent)
{
    vector<VERTEX>::iterator it;

    NEC_PQ* child_pq;
    NEC_PQ* parent_pq;

    vector<NECNode*> visited_nec_neighbors;
    VERTEX data_node_id;
    bool matched;

    /*
    cout << "CR - Data node candidates " ;
    for (int i = 0; i < data_node_candidates->size(); i += 1) {
        cout << data_node_candidates->at(i) << " ";
    }
    cout << endl;
    */
    it = data_node_candidates->begin();
    for(; it != data_node_candidates->end(); ++it) {
        data_node_id = *(it);
        if (marked_vertex->count(data_node_id) > 0 ||
            DegreeFilter(dg, data_node_id, nec_root) == true)
            //NeighborLabelFilter(dg, data_node_id, nec_root))
            continue;

        marked_vertex->insert(data_node_id);
        matched = true;

        // Need to sort children
        child_pq = GetNECEvalOrder(dg, nec_root, data_node_id, CHILD);
        parent_pq = GetNECEvalOrder(dg, nec_root, data_node_id, PARENT);

        // Recursively generate CR for child NEC, by interleaving parent and
        // child NEC nodes.
        while (child_pq->size() > 0 || parent_pq->size() > 0) {
            if(ExploreSubRegions(dg, marked_vertex, cr, parent, data_node_id, 
                                 &visited_nec_neighbors, child_pq) == false) {
                matched = false;
                break;
            }
            if (ExploreSubRegions(dg, marked_vertex, cr, parent, data_node_id, 
                                  &visited_nec_neighbors, parent_pq) == false) {
                matched = false;
                break;
            }
        }

        marked_vertex->erase(data_node_id);
        if (matched == false) {
            //cout << "Matched False :" << parent << " , " << nec_root << endl;
            continue;
        }

        // Insert into candidate region 
        InsertIntoCandidateRegion(cr, parent, nec_root, data_node_id);
    }

    //cout << "Query:" << nec_root << "_Parent:" << parent << " ";
    //cout << cr->count(nec_root) << endl;

    if (cr->count(nec_root) == 0) {
        return false;
    }

    if (cr->find(nec_root)->second->find(parent)->second->size() 
                < nec_root->members->size())
    {
        cr->find(nec_root)->second->find(parent)->second->clear();
        return false;
    }
    return true;
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

// TODO : New scoring that sums all NEC nodes not only leaf? Leaf 
// may have few candidate but intermediate NEC may have many.
ULONG MatchingOrderCalcScore(vector<NECNode*>* path, ULONG nontree_edges, CandidateRegions* cr) {
    unordered_map<VERTEX, vector<VERTEX>*>* cr_nodes;
    unordered_map<VERTEX, vector<VERTEX>*>::iterator it;
    NECNode* leaf = path->back();
    ULONG NEC_member_count = path->back()->members->size();
    ULONG score = 0;

    if (cr->count(leaf) == 0) {
        cout << "Invalid candidate region" << endl;
    }
    cr_nodes = cr->at(leaf);
    it = cr_nodes->begin();
    if ( NEC_member_count > 1) {

        while (it != cr_nodes->end()) {
            score += xChoseY(it->second->size(), NEC_member_count);
            ++ it;
        }
    } else {
        while (it != cr_nodes->end()) {
            score += it->second->size();
            ++it;
        }
        score = score / (nontree_edges + 1);
    }
    return score;
}

vector<NECNode*>* GenerateNECOrder(MatchingOrderPq* path_pq) {
    unordered_set<NECNode*> seen_nec;
    vector<NECNode*>* nec_queue;
    MatchingOrderPair* pair;
    NECNode* curr_nec;
    vector<NECNode*>* evaluation_queue = new vector<NECNode*>;

    while (!path_pq->empty()) {
        pair = path_pq->top();
        path_pq->pop();
        nec_queue = pair->q;
        
        for (int index = nec_queue->size() - 1; index >= 0; index -= 1) {
            curr_nec = nec_queue->at(index);
            if (seen_nec.count(curr_nec) > 0)
                continue;

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

MatchingOrderPq* MatchingOrder(NECNode* nec_root, CandidateRegions* cr) {
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
            pq_elem->score = MatchingOrderCalcScore(path, nontree_edges, cr);
            pq->push(pq_elem);

            /*
            cout << "****** copied queue size: " << path->size() << endl;
            for (int i = 0; i < path->size(); i += 1) {
                cout << path->at(i) << endl;
            }
            cout << "printed copied queue" << endl;
            */
        }
    }

    return pq;
}

CandidateRegions* AllocCandidateRegions(DataGraph& dg, NECNode* nec_root, LABEL root_label) {
    CandidateRegions* cr = new CandidateRegions;
    unordered_set<VERTEX>* visited_vertices = new unordered_set<VERTEX>;
    VERTEX root_node_id;
    vector<VERTEX> candidate_data_vertices;
    VERTEX* data_node_iter;

    data_node_iter = dg.vfs->GetVertexIterator(root_label);
    for(int data_node_index = 0; 
        data_node_index < dg.vfs->NumVertices(root_label); 
        data_node_index += 1)
    {
        // TODO : FILTER data nodes based on NEC characteristics.
        root_node_id = *(data_node_iter + data_node_index);
        candidate_data_vertices.clear();
        candidate_data_vertices.push_back(root_node_id);

        if(ExploreCandidateRegions(dg, nec_root, &candidate_data_vertices, 
                                   visited_vertices, cr, -1) == false)
        {
            continue;
        }

        /*
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
        }*/
        // Iterate and print me
    }
    return cr;
}

bool SubGraphChildrenReferences(vector<VERTEX>* children, 
                                VERTEX data_node,
                                DataGraph& dg) {

    VERTEX child;
    for (int i = 0; i < children->size(); i += 1) {
        child = children->at(i);
        if (!dg.efs->HasEdge(data_node, child)) {
            //cout << "no edge from " << data_node << " to " << child << endl;
            return false;
        }
    }

    return true;
}

vector<VERTEX>* SubGraphDataVertexUnion(vector<vector<VERTEX>*>* collections,
                                        unordered_set<VERTEX>* marked_data_vertex,
                                        vector<VERTEX>* matched_child_dnodes,
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
                SubGraphChildrenReferences(matched_child_dnodes, data_node, dg))
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

vector<VERTEX>* SubGraphDataVertexIntersection(vector<vector<VERTEX>*>* collections,
                                               unordered_set<VERTEX>* marked_data_vertex,
                                               vector<VERTEX>* matched_child_dnodes,
                                               DataGraph& dg) 
{
    if (collections->size() == 0)
        return NULL;

    vector<VERTEX>* seed = collections->at(0);
    VERTEX matched_node;
    unordered_set<VERTEX>* result = new unordered_set<VERTEX>();
    unordered_set<VERTEX>::iterator it;

    for (int i = 0; i < seed->size(); i += 1) {
        matched_node = seed->at(i);
        if (marked_data_vertex->count(matched_node) == 0 &&
            SubGraphChildrenReferences(matched_child_dnodes, matched_node, dg))
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
                it++;
            } else {
                it = result->erase(it);
            }
        }
    }
    return NULL;
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



bool SubGraphCountCRCmp(vector<VERTEX>* v1, vector<VERTEX>* v2) {
    return (v1->size() < v2->size());
}

// Get candidate data vertex for a NEC node
vector<VERTEX>* SubGraphCountGetDataVertex(NECNode* nec, 
                                unordered_map<QueryNode*, VERTEX>* mapping,
                                unordered_set<NECNode*>* matched_nec,
                                unordered_set<VERTEX>* marked_data_vertex,
                                CandidateRegions* cr,
                                DataGraph& dg)
{
    NECNode* neighbor;
    QueryNode* query_node;
    VERTEX mapped_data_node;
    unordered_map<VERTEX, vector<VERTEX>*>* nec_cr_map;
    unordered_map<VERTEX, vector<VERTEX>*>::iterator cr_iter;
    vector<vector<VERTEX>*> valid_candidate_regions;
    vector<VERTEX>* matched_child_dnodes;

    // NEC doesn't exist. Means there was no valid CR.
    if (cr->count(nec) == 0) {
        return NULL;
    }

    nec_cr_map = cr->at(nec);
    matched_child_dnodes = SubGraphGetMatchedChildren(nec, matched_nec, mapping);

    /*
    cout << "MATCHED CHILD QUERY NODES " ;
    for (int i = 0; i < matched_child_dnodes->size(); i += 1) {
        cout << matched_child_dnodes->at(i) << " - ";
    }
    cout << endl;
    */
    // Find candidate data vertex based on pre-matched query vertex
    for (int i = 0; i < nec->parents->size(); i += 1) {
        neighbor = nec->parents->at(i);
        if (matched_nec->count(neighbor) == 0)
            continue;

        for (int member_index = 0;
             member_index < neighbor->members->size();
             member_index += 1)
        {
            query_node = neighbor->members->at(member_index);
            // If NEC has been matched, all of it's query vertex must be
            // matched as well.
            mapped_data_node = mapping->at(query_node);
            if (nec_cr_map->count(mapped_data_node) > 0) {
                valid_candidate_regions
                    .push_back(nec_cr_map->at(mapped_data_node));
            }
        }
    }

    // If no candidate data vertex were found, means no neighbor NEC are matched
    // so take CR(u,*) for candidate data vertex
    if (valid_candidate_regions.size() == 0) {
        cr_iter = nec_cr_map->begin();
        while (cr_iter != nec_cr_map->end()) {
            valid_candidate_regions.push_back(cr_iter->second);
            cr_iter ++;
        }
        return SubGraphDataVertexUnion(&valid_candidate_regions, 
                                       marked_data_vertex, 
                                       matched_child_dnodes,
                                       dg);
    } else {
        // Sort CR vector by size for join order.
        sort (valid_candidate_regions.begin(),
              valid_candidate_regions.end(),
              SubGraphCountCRCmp);
        return SubGraphDataVertexIntersection(&valid_candidate_regions, 
                                              marked_data_vertex,
                                              matched_child_dnodes,
                                              dg);
    }
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

int SubGraphCount(vector<NECNode*>* matching_order, int matching_index,
                   unordered_map<QueryNode*, VERTEX>* mapping, 
                   unordered_set<NECNode*>* matched_nec,
                   unordered_set<VERTEX>* marked_data_vertex,
                   CandidateRegions* cr,
                   DataGraph& dg) 
{
    NECNode* curr_nec;
    vector<VERTEX>* candidate_data_vertex;
    vector<VERTEX>* combination;
    bool matched;
    int count = 0;
    curr_nec = matching_order->at(matching_index);
    candidate_data_vertex = SubGraphCountGetDataVertex(curr_nec, mapping,
                                                       matched_nec,
                                                       marked_data_vertex,
                                                       cr,
                                                       dg);

    /*
    cout << "NEC:" << curr_nec << " CANIDATE SIZE" << candidate_data_vertex->size() <<  endl;
    cout << "DATA VERTEX CANDIDATE:";
    for (int i = 0; i < candidate_data_vertex->size(); i += 1) {
        cout << candidate_data_vertex->at(i) << " - ";
    }
    cout << endl;
    */
    XChoseY combinations = XChoseY(candidate_data_vertex, curr_nec->members->size());
    combination = combinations.getNext();

    while (combination != NULL) {
        /*
        cout << "Combination :";
        for (int i = 0; i < combination->size(); i += 1)
            cout << combination->at(i) << " - ";
        cout << endl;
        */

        matched = true;
        for (int i = 0; i < combination->size(); i += 1) {
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
            count += 1;
        } else {
            count += SubGraphCount(matching_order, matching_index + 1, 
                                   mapping, matched_nec, 
                                   marked_data_vertex, cr, dg);
        }
        SubGraphRestore(mapping, combination, curr_nec, matched_nec, 
                        marked_data_vertex);

        combination = combinations.getNext();
    }
    return count;
}

void TurboIso(DataGraph& dg, query_node_map* qg) {
    QueryNode* start_vertex;
    LABEL root_label;
    NECNode* nec_root;
    CandidateRegions* cr;
    MatchingOrderPq* matching_order_pq;
    vector<NECNode*>* matching_queue;
    unordered_set<NECNode*> matched_NEC;
    unordered_map<QueryNode*, VERTEX> query_node_mapping;
    unordered_set<VERTEX> marked_data_vertex;

    cout << "start vertex" << endl;
    // Get starting query vertex for NEC tree construction
    start_vertex = GetStartingQueryVertex(dg, qg);


    cout << "nec tree" << endl;
    // Construct NEC tree and populate CR
    root_label = start_vertex->label;
    nec_root = RewriteToNECTree(qg, start_vertex);
    cout << "cr" << endl;
    cr = AllocCandidateRegions(dg, nec_root, root_label);
    //SortCr(cr);
    //
    cout << "match order" << endl;
    // Get Matching Order for NEC
    matching_order_pq = MatchingOrder(nec_root, cr);
    matching_queue = GenerateNECOrder(matching_order_pq);
    
    cout << "count subgraph" << endl;
    // Count SubGraphs
    matched_NEC.clear();
    query_node_mapping.clear();
    int answer = SubGraphCount(matching_queue, 0, &query_node_mapping, &matched_NEC,
                  &marked_data_vertex, cr, dg);
    cout << "COUNT: " << answer << endl;
}

int main() {
    DataGraph dg;
    query_node_map* qg;

    VertexFs vfs("dblp_vfs");
    EdgeFs efs("dblp_efs.txt", "dblp_efs_rev.txt", &vfs);
    dg.vfs = &vfs;
    dg.efs = &efs;

    qg = ReadQueryGraphFromFile("query_graph.txt");
    /*
    cout << "print query graph" << endl;
    PrintQueryGraph(qg);
    */

    clock_t begin = clock();
    TurboIso(dg, qg);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "secs:" << elapsed_secs << endl;

}
#endif
