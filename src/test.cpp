#include <iostream>
using namespace std;
#include "combinators.cpp"
#include "assert.h"



void checkCase(vector<VERTEX>* combination, VERTEX val1, VERTEX val2) {

    cout << "Combination:" ;
    cout << (VERTEX)combination->at(0) << " - ";
    cout << (VERTEX)combination->at(1) << endl;
    
    assert(val1 == (VERTEX)combination->at(0));
    assert(val2 == (VERTEX)combination->at(1));
    
}
void check3Case(vector<VERTEX>* combination, VERTEX val1, VERTEX val2, VERTEX val3) {

    cout << "Combination:" ;
    cout << (VERTEX)combination->at(0) << " - ";
    cout << (VERTEX)combination->at(1) << " - ";
    cout << (VERTEX)combination->at(2) << endl;
    assert(val1 == (VERTEX)combination->at(0));
    assert(val2 == (VERTEX)combination->at(1));
    assert(val3 == (VERTEX)combination->at(2));
}

void testXChoseY() {
    vector<VERTEX> vals = vector<VERTEX>(4, 0);
    vector<VERTEX>* result;
    for (int i = 0; i < 4; i += 1) {
        vals[i] = i;
    }

    XChoseY obj = XChoseY((vector<VERTEX>*) &vals, 2);
    // 6 cases
    checkCase(obj.getNext(), 0, 1);
    checkCase(obj.getNext(), 0, 2);
    checkCase(obj.getNext(), 0, 3);
    checkCase(obj.getNext(), 1, 2);
    checkCase(obj.getNext(), 1, 3);
    checkCase(obj.getNext(), 2, 3);

    assert(obj.getNext() == NULL);
}

void testManyChoseOne() {
    vector<VERTEX> bucket1 = vector<VERTEX>(1,0);
    bucket1[0] = 1;
    vector<VERTEX> bucket2 = vector<VERTEX>(3,0);
    bucket2[0] = 2;
    bucket2[1] = 4;
    bucket2[2] = 3;
    vector<VERTEX> bucket3 = vector<VERTEX>(2,0);
    bucket3[0] = 5;
    bucket3[1] = 6;

    ManyChoseOne obj;
    obj.addBucket((vector<VERTEX>*) &bucket1);
    obj.addBucket((vector<VERTEX>*) &bucket2);
    obj.addBucket((vector<VERTEX>*) &bucket3);

    check3Case(obj.getNext(), 1, 2, 5);
    check3Case(obj.getNext(), 1, 2, 6);
    check3Case(obj.getNext(), 1, 4, 5);
    check3Case(obj.getNext(), 1, 4, 6);
    check3Case(obj.getNext(), 1, 3, 5);
    check3Case(obj.getNext(), 1, 3, 6);
}

int main() {
    cout << "Test" << endl;
    testXChoseY();
    testManyChoseOne();
    cout << "Done" << endl;
}
