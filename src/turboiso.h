#ifndef TURBOISO_H
#define TURBOISO_H

struct RankingPair {
    QueryNode* v;
    double score;
};

/*
 *  Sorts for min-pq
 */
struct CompareRankingPair {
    bool operator()(const RankingPair* p1, const RankingPair* p2) {
        return p1->score > p2->score;
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
#endif
