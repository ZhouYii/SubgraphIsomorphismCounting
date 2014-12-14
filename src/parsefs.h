#ifndef PARSEFS_H_
#define PARSEFS_H_

#define UINT unsigned int
#define VERTEX unsigned int
#define LABEL unsigned int
#define ULONG unsigned long

#include "graph.cpp"
#include <iostream>
#include <queue>
#include <utility>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
using namespace std;
typedef std::unordered_map<UINT, QueryNode*> query_node_map;


class VertexFs 
{
    private :
        UINT* data;
        UINT size;

        UINT* label_names_offset;
        UINT* label_end_offset;
        UINT* vertex_offset;

    public :
        VertexFs(const char* file_path);

        UINT NumVertices();
        UINT NumVertices(const UINT label);
        UINT NumLabels();

        UINT* GetVertexIterator();
        UINT* GetVertexIterator(const UINT label);
        LABEL GetLabel(const UINT vertex);

        void WriteToFile(const char* file_path);

};

class EdgeFs
{
    private :
        UINT size;
        UINT reverse_size;
        UINT num_vertices;
        // Offsets for r_boundaries of end_offsets reflect the fact
        // that each (label,offset) pair in adj_vertex_labels is 2 UINTS
        UINT* data;
        UINT* reverse_data;

        // Start of Block 2: (Label,Offset) Pair block
        UINT* adjacent_vertex_labels;
        // Start of Block 1 : VertexId block
        UINT* adjacent_vertex_ids;
        UINT* rev_adjacent_vertex_labels;
        UINT* rev_adjacent_vertex_ids;

        VertexFs* vfs;

        unordered_map<VERTEX,unordered_map<LABEL, vector<VERTEX>*>*>* 
                        GenerateReverseIndex(UINT& size);
        void ConstructReverseEdges();

    public :

        EdgeFs(const char* file_path, const char* rev_file_path, VertexFs* vfs);

        // Mutators
        void FlipEdges();

        // Return iterator to vertex_id block
        UINT* OutgoingEdgeIterator(const UINT vertex_id);
        UINT* IncomingEdgeIterator(const UINT vertex_id);
        // Need to handle NULL
        UINT* OutgoingEdgeIterator(const UINT vertex_id, const LABEL label);
        UINT* OutgoingEdgeIteratorEnd(const UINT vertex_id, const LABEL label);
        UINT* IncomingEdgeIterator(const UINT vertex_id, const LABEL label);

        UINT* OutgoingLabelIterator(const UINT vertex_id);
        UINT* IncomingLabelIterator(const UINT vertex_id);

        vector<VERTEX>* GetOutdegreeVertices(VERTEX src, LABEL l);
        vector<VERTEX>* GetIndegreeVertices(VERTEX src, LABEL l);
        UINT GetInDegree(const UINT vertex_id);
        UINT GetOutDegree(const UINT vertex_id);
        UINT GetInDegree(const UINT vertex_id, const LABEL l);
        UINT GetOutDegree(const UINT vertex_id, const LABEL l);

        bool HasEdge(const VERTEX from, const VERTEX to);


        // Human readable,editable
        void PrintGraphEdges(const char* file_path);
        void WriteToFile(const char* file_path);
};

struct DataGraph {
    VertexFs* vfs;
    EdgeFs* efs;
};

/*
 * Used for PQ comparator
 */
struct CompareLabels {
    bool operator()(const pair<LABEL, vector<VERTEX>*>* it,
                    const pair<LABEL, vector<VERTEX>*>* other) {
        return it->first > other->first;
    }
};

#endif
