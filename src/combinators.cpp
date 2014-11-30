#ifndef combinator_cpp
#define combinator_cpp

#include "defines.h"
#include "combinators.h"

XChoseY::XChoseY(vector<NECNode*>* bucket, int num_to_select) {
    this->bucket = bucket;
    this->chose_num = num_to_select;
    this->iterators = new vector<int>(num_to_select, 0);
    for (int i = 0; i < num_to_select; i += 1) {
        this->iterators->at(i) += i;
    }
}

XChoseY::~XChoseY() {
    if (iterators != NULL)
        delete iterators;
}

vector<NECNode*>* XChoseY::getNext() {
    vector<NECNode*>* result;
    int index;
    int last_advanced_iterator;

    // Done case.
    if (iterators == NULL) {
        return NULL;
    }

    result = new vector<NECNode*>(chose_num, 0);
    for (int iter_num = 0; iter_num < chose_num; iter_num += 1) {
        index = this->iterators->at(iter_num);
        result->at(iter_num) = this->bucket->at(index);
    }

    last_advanced_iterator = advanceOneIterator(chose_num - 1);
    if (last_advanced_iterator == -1) {
        // Done. Free mem.
        delete iterators;
        iterators = NULL;
    }
    return result;
}

// Index is which iterator index is being incremented.
int XChoseY::advanceOneIterator(int index) {
    int curr_index;
    int dst_from_end;
    int max_index;
    int modified_iterator_number;
    if (iterators == NULL)
        return -1;

    modified_iterator_number = index;
    curr_index = iterators->at(index);
    dst_from_end = iterators->size() - index - 1;
    max_index = bucket->size() - dst_from_end - 1;

    if (curr_index == max_index) {
        // If 0th iterator reaches max index, no more combinations exist
        if (index > 0) {
            modified_iterator_number = advanceOneIterator(index - 1);
            iterators->at(index) = iterators->at(index - 1) + 1;
        } else {
            modified_iterator_number = -1;
        }
    } else {
        iterators->at(index) += 1;
    }

    return modified_iterator_number;
}

ManyChoseOne::ManyChoseOne() {
    this->buckets = new vector<vector<NECNode*>*>;
    this->iterators = new vector<int>;
    this->has_remaining = true;
}

ManyChoseOne::~ManyChoseOne() {
    delete buckets;
    delete iterators;
}

void ManyChoseOne::addBucket(vector<NECNode*>* bucket) {
    buckets->push_back(bucket);
    iterators->push_back(0);
}

bool ManyChoseOne::incrementBucketIterator(int bucket_num) {
    if (bucket_num < 0) {
        this->has_remaining = false;
        return true;
    }
    bool result = false;
    int max_index = buckets->at(bucket_num)->size() - 1;
    int curr_index = iterators->at(bucket_num);

    if (curr_index < max_index) {
        iterators->at(bucket_num) += 1;
    } else {
        iterators->at(bucket_num) = 0;
        incrementBucketIterator(bucket_num - 1);
        result = true;
    }

    return result;
}

NECNode* ManyChoseOne::getCurrentElemFromBucket(int bucket_index) {
    int index = iterators->at(bucket_index);
    vector<NECNode*>* bucket = buckets->at(bucket_index);
    return bucket->at(index);
}

vector<NECNode*>* ManyChoseOne::getNext() {
    unordered_set<NECNode*> seen;
    vector<NECNode*>* bucket;
    int bucket_index;
    int bucket_num;
    bool done = false;
    bool bucket_index_reset;
    NECNode* curr_nec;
    vector<NECNode*>* result;

    if (!this->has_remaining) {
        return NULL;
    }

    result = new vector<NECNode*>(buckets->size(), 0);
    while (!done && this->has_remaining) {
        for (bucket_num = 0; bucket_num < buckets->size(); bucket_num += 1) {
            bucket_index_reset = false;
            curr_nec = getCurrentElemFromBucket(bucket_num);
            if (seen.count(curr_nec) > 0) {
                do {
                    bucket_index_reset = incrementBucketIterator(bucket_num);
                    // Overwrite accumulated results because we reset the current
                    // iterator as well as previous bucket's iterators
                    if (bucket_index_reset) {
                        seen.clear();
                        break;
                    }

                    curr_nec = getCurrentElemFromBucket(bucket_num);
                } while (seen.count(curr_nec) > 0);
            }
            seen.insert(curr_nec);
            result->at(bucket_num) = curr_nec;
        }

        // If for loop completed, result is valid. Otherwise, restart.
        if (bucket_num >= buckets->size()) {
            incrementBucketIterator(buckets->size() - 1);
            done = true;
        }
     }
    return result;
}

#endif
