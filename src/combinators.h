#ifndef combinator_h
#define combinator_h

class XChoseY {
    private :
        // Should not free this.
        vector<NECNode*>* bucket;
        vector<int>* iterators;
        int chose_num;

        int advanceOneIterator(int index);

    public :
        XChoseY(vector<NECNode*>* bucket, int num_to_select);
        ~XChoseY();
        vector<NECNode*>* getNext();
};

// Many buckets, chose one from each. Unique.
class ManyChoseOne {
    private :
        bool has_remaining;
        vector<vector<NECNode*>*>* buckets;
        vector<int>* iterators;

        bool incrementBucketIterator(int bucket_index);
        NECNode* getCurrentElemFromBucket(int bucket_index);

    public :
        ManyChoseOne();
        ~ManyChoseOne();
        vector<NECNode*>* getNext();
        void addBucket(vector<NECNode*>* bucket);
};

#endif
