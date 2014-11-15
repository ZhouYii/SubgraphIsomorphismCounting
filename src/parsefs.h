#ifndef PARSEFS_H_
#define PARSEFS_H_

#include "defines.h"
#include "graph.cpp"
#include <iostream>
#include <queue>
#include <utility>
#include <stdio.h>
using namespace std;

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
        // Offsets for r_boundaries of adjacent_vertex_labels reflect the fact
        // that each (label,offset) pair is 2 UINTS
        UINT* data;
        UINT* reverse_data;
        UINT* adjacent_vertex_labels;
        UINT* adjacent_vertex_ids;

        VertexFs* vfs;

        unordered_map<VERTEX,unordered_map<LABEL, vector<VERTEX>*>*>* 
                        GenerateReverseIndex(UINT& size);
        void ConstructReverseEdges();

    public :
        EdgeFs(const char* file_path, VertexFs* vfs);

        // Mutators
        void FlipEdges();

        UINT* GetEdgeIterator(const UINT vertex_id);
        UINT GetOutDegree(const UINT vertex_id);

        // Human readable,editable
        void PrintGraphEdges(const char* file_path);
        void WriteToFile(const char* file_path);
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
