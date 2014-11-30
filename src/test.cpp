#include <iostream>
using namespace std;
#include "combinators.cpp"
#include "assert.h"



void checkCase(vector<NECNode*>* combination, long val1, long val2) {
    assert(val1 == (long)combination->at(0));
    assert(val2 == (long)combination->at(1));
}
void check3Case(vector<NECNode*>* combination, long val1, long val2, long val3) {

    cout << "Combination:" ;
    cout << (long)combination->at(0) << " - ";
    cout << (long)combination->at(1) << " - ";
    cout << (long)combination->at(2) << endl;
    assert(val1 == (long)combination->at(0));
    assert(val2 == (long)combination->at(1));
    assert(val3 == (long)combination->at(2));
}

void testXChoseY() {
    vector<long> vals = vector<long>(4, 0);
    vector<NECNode*>* result;
    for (int i = 0; i < 4; i += 1) {
        vals[i] = i;
    }

    XChoseY obj = XChoseY((vector<NECNode*>*) &vals, 2);
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
    vector<long> bucket1 = vector<long>(1,0);
    bucket1[0] = 1;
    vector<long> bucket2 = vector<long>(3,0);
    bucket2[0] = 2;
    bucket2[1] = 4;
    bucket2[2] = 3;
    vector<long> bucket3 = vector<long>(2,0);
    bucket3[0] = 5;
    bucket3[1] = 6;

    ManyChoseOne obj;
    obj.addBucket((vector<NECNode*>*) &bucket1);
    obj.addBucket((vector<NECNode*>*) &bucket2);
    obj.addBucket((vector<NECNode*>*) &bucket3);

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
