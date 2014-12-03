#ifndef NODE_DEFINITIONS_CPP
#define NODE_DEFINITIONS_CPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#define UINT unsigned int

using namespace std;

struct QueryNode {
    UINT label;
    vector<QueryNode*>* parents;
    vector<QueryNode*>* children;
};

struct NECNode {
    // Members implicitly represent NEC attributes
    vector<QueryNode*>* members;

    // Storage for directional edges
    vector<NECNode*>* parents;
    vector<NECNode*>* children;
};

typedef std::unordered_map<UINT, QueryNode*> query_node_map;

struct QueryGraph {
    query_node_map* query_nodes;
    unordered_map<QueryNode*, UINT>* initialized_query_nodes;
};

NECNode* NECNodeInit() {
    NECNode* nec_node = new NECNode;
    nec_node->members = new vector<QueryNode*>;
    nec_node->children = new vector<NECNode*>;
    nec_node->parents = new vector<NECNode*>;
    return nec_node;
}

vector<string>& SplitString(const string &s, char delim, vector<string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void PrintQueryGraph(query_node_map* qm) {
    query_node_map::iterator it = qm->begin();
    while (it != qm->end()) {
        cout << "ID:" << it->second;
        cout << " Children:";
        for (int i = 0; i < it->second->children->size(); i += 1) {
            cout << " - " << it->second->children->at(i);
        }
        cout << " Parents:";
        for (int i = 0; i < it->second->parents->size(); i += 1) {
            cout << " - " << it->second->parents->at(i);
        }
        cout << endl;
        ++it;
    }
}

QueryGraph* ReadQueryGraphFromFile(const char* file_path) {
    QueryGraph* qg = new QueryGraph;
    query_node_map* map = new query_node_map;
    unordered_map<QueryNode*, UINT>* initialized = new unordered_map<QueryNode*, UINT>;
    QueryNode* vertex;
    string line;
    vector<string> str_buf;

    qg->query_nodes = map;
    qg->initialized_query_nodes = initialized;
    ifstream file(file_path);
    // Parse node id and labels
    while (getline(file,line)) {
        if (line[0] == '#')
            continue;
        if (line.length() == 0)
            break;

        SplitString(line, '\t', str_buf);
        vertex = new QueryNode;
        vertex->label = stoul(str_buf[1]);
        vertex->parents = new vector<QueryNode*>;
        vertex->children = new vector<QueryNode*>;
        map->insert(query_node_map::value_type(stoul(str_buf[0]), vertex));
        str_buf.clear();

        cout << "QueryNode Mapping :"<< stoul(str_buf[0]) << "-";
        cout << vertex << endl;
    }
    // Parse adjacency list
    while (getline(file,line)) {
        if (line[0] == '#') {
            continue;
        }

        if (line.length() == 0)
            break;

        SplitString(line, '\t', str_buf);
        vertex = map->find(stoul(str_buf[0]))->second;
        for(int i = 2; i < str_buf.size(); i += 1) {
            vertex->children->push_back(map->find(stoul(str_buf[i]))->second);
            map->find(stoul(str_buf[i]))->second->parents->push_back(vertex);
        }
        str_buf.clear();
    }
    // Parse initialized query vertices.
    while (getline(file,line)) {
        if (line[0] == '#') {
            continue;
        }

        if (line.length() == 0)
            break;

        SplitString(line, '\t', str_buf);

        UINT qnode_id = stoul(str_buf[0]);
        UINT dnode_id = stoul(str_buf[1]);
        if (map->count(qnode_id) > 0) {
            QueryNode* qnode = map->at(qnode_id);
            initialized->insert(make_pair<QueryNode*, UINT>(qnode, dnode_id));
        }
        str_buf.clear();
    }

    cout << "initialized size: " << initialized->size() << endl;
    unordered_map<QueryNode*,UINT>::iterator it = initialized->begin();
    while (it != initialized->end()) {
        cout << it->first << " " << it->second << endl;
        ++it;
    }


    return qg;
}
/*
int main() {
    ReadQueryGraphFromFile("query_graph.txt");
}
*/
#endif
