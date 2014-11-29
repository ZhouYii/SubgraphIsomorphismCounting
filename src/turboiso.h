#ifndef TURBOISO_H
#define TURBOISO_H

struct MatchingOrderPair {
    queue<NECNode*>* q;
    ULONG score;
};
struct RankingPair {
    QueryNode* v;
    double score;
};

/*
 *  Sorts for min-pq. TODO : verify min pq not max pq
 */
struct CompareRankingPair {
    bool operator()(const RankingPair* p1, const RankingPair* p2) {
        return p1->score > p2->score;
    }
};

/*
 * Min PQ for matching order.
 * */
struct CompareMatchingOrder {
    bool operator()(const MatchingOrderPair* m1, const MatchingOrderPair* m2) {
        return m1->score > m2->score;
    }
};

/*
 * Sorting for NEC eval order
 */
struct CompareNEC {
    bool operator()(const pair<NECNode*, UINT>* p1, pair<NECNode*, UINT>* p2) {
        return p1->second < p2->second;
    }
};

/*
 * Used when BFS-ing through query graph.
 * A node's parent might not be reachable from BFS on the node's chidren alone.
 * When children BFS queue runs out, BFS through parent's BFS queue.
 */
struct AncestorEquivalentClass {
    QueryNode* v;
    QueryNode* child;
};

struct CRNode {
    // Encapsulates the 'v' prior node
    vector<VERTEX> candidate_data_nodes;
    // Ptrs to CRNode of adjacent NEC vertex, with current CRNode assignment as
    // a prior (2nd level map has domain = current CR's data vertexx)
    unordered_map<NECNode*, unordered_map<VERTEX, CRNode*> > children;
    unordered_map<NECNode*, unordered_map<VERTEX, CRNode*> > parents;
};

struct CRTree {
    CRNode* root;
    unordered_map<NECNode*, unordered_set<CRNode*> > candidate_regions;
};

typedef priority_queue<MatchingOrderPair*,
                       vector<MatchingOrderPair*>,
                       CompareMatchingOrder> MatchingOrderPq;
#endif
