#ifndef combinator_h
#define combinator_h

class XChoseY {
    private :
        // Should not free this.
        vector<VERTEX>* bucket;
        vector<int>* iterators;
        int chose_num;

        int advanceOneIterator(int index);

    public :
        XChoseY(vector<VERTEX>* bucket, int num_to_select);
        ~XChoseY();
        vector<VERTEX>* getNext();
};

// Many buckets, chose one from each. Unique.
class ManyChoseOne {
    private :
        bool has_remaining;
        vector<vector<VERTEX>*>* buckets;
        vector<int>* iterators;

        bool incrementBucketIterator(int bucket_index);
        VERTEX getCurrentElemFromBucket(int bucket_index);

    public :
        ManyChoseOne();
        ~ManyChoseOne();
        vector<VERTEX>* getNext();
        void addBucket(vector<VERTEX>* bucket);
};

#endif
