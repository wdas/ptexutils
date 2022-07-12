/*
* Copyright Disney Enterprises, Inc.  All rights reserved.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License
* and the following modification to it: Section 6 Trademarks.
* deleted and replaced with:
*
* 6. Trademarks. This License does not grant permission to use the
* trade names, trademarks, service marks, or product names of the
* Licensor and its affiliates, except as required for reproducing
* the content of the NOTICE file.
*
* You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef KdTree_h
#define KdTree_h

/* balanced kdtree
   Brent Burley, Mar 2006

   The tree is left-balanced (i.e. the left subtree is always
   complete).  The nodes are stored in traversal order ("pre-order")
   to optimize memory lookahead and cache-coherence during traversal.

   This means that:
   1) The root node is node 0, depth 0.
   2) Any given subtree is stored in a contiguous memory block.
   3) The left child of node n is node n+1.
   4) The right child of node n is node n+(size of left subtree).
   5) The size of the left subtree of any node can be easily
      determined based on the node's overall subtree size (left+right+1).
      This can be propagated down during traversal.
*/

#include <string.h>
#include <vector>
#include <float.h>
#include <algorithm>
#include <cassert>

template <int k> class BBox
{
 public:
    float min[k];
    float max[k];

    BBox() { clear(); }
    BBox(const float p[k]) { set(p); }

    void set(const float p[k])
    {
        for (int i = 0; i < k; i++) {
            min[i] = max[i] = p[i];
        }
    }

    void clear()
    {
        for (int i = 0; i < k; i++) {
            min[i] = FLT_MAX;
            max[i] = FLT_MIN;
        }
    }


    void grow(const float p[k])
    {
        for (int i = 0; i < k; i++) {
            if (p[i] < min[i]) min[i] = p[i];
            if (p[i] > max[i]) max[i] = p[i];
        }
    }

    void grow(float R)
    {
        for (int i = 0; i < k; i++) {
            min[i] -= R;
            max[i] += R;
        }
    }

    void grow(const BBox& b)
    {
        for (int i = 0; i < k; i++) {
            if (b.min[i] < min[i]) min[i] = b.min[i];
            if (b.max[i] > max[i]) max[i] = b.max[i];
        }
    }

    bool intersects(const BBox& b) const
    {
        for (int i = 0; i < k; i++) {
            if (b.max[i] < min[i] || b.min[i] > max[i]) return 0;
        }
        return 1;
    }

    bool inside(const float p[k]) const
    {
        for (int i = 0; i < k; i++) {
            if (p[i] < min[i] || p[i] > max[i]) return 0;
        }
        return 1;
    }
};


// MAKE-heap on a binary (array based) heap
inline float buildHeap(std::vector<int>& result,std::vector<float>& distance_squared)
{
    assert(result.size()==distance_squared.size());
    int heap_size=result.size();
    int max_non_leaf_index=result.size()/2-1; // 0 to max_non_leaf_index is indices of parents

    // go from bottom of tree scanning right to left in each height of the tree
    for (int subtreeParent=max_non_leaf_index;subtreeParent>=0;subtreeParent--){ 
        int current_parent=subtreeParent;
        while (current_parent<=max_non_leaf_index){
            int left_index=2*current_parent+1;int right_index=2*current_parent+2;
            // find largest element
            int largest_index=current_parent;
            if (left_index<heap_size && distance_squared[left_index]>distance_squared[largest_index])
                largest_index=left_index;
            if (right_index<heap_size && distance_squared[right_index]>distance_squared[largest_index])
                largest_index=right_index;
            // subtree parented at current_parent satisfies heap property because parent has largest
            if (largest_index==current_parent) break;
            // pull large value up and descend to tree to see if new that subtree is invalid
            std::swap(result[largest_index],result[current_parent]);
            std::swap(distance_squared[largest_index],distance_squared[current_parent]);
            current_parent=largest_index;
        }
    }
    return distance_squared[0]; // return max distance
}

// Inserts smaller element into heap (does not check so caller must)
inline float insertToHeap(std::vector<int>& result,std::vector<float>& distance_squared,int new_id,float new_distance_squared)
{
    assert(result.size()>0 && distance_squared.size()==result.size());
    assert(new_distance_squared<distance_squared[0]);
    int current_parent=0;
    int heap_size=result.size();
    for(;;){
        int left=2*current_parent+1,right=2*current_parent+2;
        int index_of_largest=current_parent;
        // find out if the current element is being replaced by descendant or by new element
        if (left>=heap_size) break;
        else if (right>=heap_size || distance_squared[left]>distance_squared[right]) index_of_largest=left;
        else index_of_largest=right;
        // new element is largest
        if (new_distance_squared>distance_squared[index_of_largest]) break;
        // pull the biggest element up and recurse to where it came from
        std::swap(result[index_of_largest],result[current_parent]);
        std::swap(distance_squared[index_of_largest],distance_squared[current_parent]);
        current_parent=index_of_largest;
    }
    // overwrite current node in tree
    distance_squared[current_parent]=new_distance_squared;
    result[current_parent]=new_id;
    return distance_squared[0]; // return max distance
}


template <int k> class KdTree
{

    struct NearestNQuery
    {
        NearestNQuery(std::vector<int>& result,std::vector<float>& distanceSquared,const float pquery_in[k],uint32_t maxPoints,float maxRadiusSquared)
            :result(result),distanceSquared(distanceSquared),maxPoints(maxPoints),maxRadiusSquared(maxRadiusSquared)
        {for(int i=0;i<k;i++) pquery[i]=pquery_in[i];}

        std::vector<int>& result;
        std::vector<float>& distanceSquared;
        float pquery[k];
        uint32_t maxPoints;
        float maxRadiusSquared;
    };

    struct NearestQuery
    {
        NearestQuery(const float pquery_in[k], float maxRadiusSquared)
            : nearest(-1), maxRadiusSquared(maxRadiusSquared)
        {for(int i=0;i<k;i++) pquery[i]=pquery_in[i];}

        float pquery[k];
        int nearest;
        float maxRadiusSquared;
    };

 public:
    KdTree();
    ~KdTree();
    void clear();
    inline int size() const { return _points.size(); }
    inline const BBox<k>& bbox() const { return _bbox; }
    inline void setBounds(BBox<k>& b) { _bbox = b; }
    inline void setSorted(bool b) { _sorted = b; }
    inline const float* point(int i) const { return _points[i].p; }
    inline int id(int i) const { return i<int(_ids.size()) ? _ids[i] : i; }
    void setPoints(float* p, int n);
    void allocatePoints(int n);
    void sort();
    void findPoints(std::vector<int>& points, const BBox<k>& bbox) const;
    float findNPoints(std::vector<int>& result,std::vector<float>& distanceSquared,
                      const float p[k],int nPoints,float maxRadius) const;
    int findNearest(const float p[k],float maxRadius) const;

    inline void setPoint(int i, float* p, bool updateBounds = true)
    {
        memcpy(&_points[i], p, sizeof(Point));
        if (updateBounds) _bbox.grow(_points[i].p);
    }

 private:
    void sortSubtree(int n, int count, int j);
    struct ComparePointsById {
        float* points;
        ComparePointsById(float* p) : points(p) {}
        bool operator() (int a, int b) { return points[a*k] < points[b*k]; }
    };
    void findPoints(std::vector<int>& result, const BBox<k>& bbox,
                    int n, int size, int j) const;
    void findNPoints(NearestNQuery& query,int n,int size,int j) const;
    void findNearest(NearestQuery& query,int n,int size,int j) const;

    static inline void ComputeSubtreeSizes(int size, int& left, int& right)
    {
        left = size/2;
        right = size-left-1;
    } 

    BBox<k> _bbox;
    struct Point { float p[k]; };
    std::vector<Point> _points;
    std::vector<int> _ids;
    bool _sorted;
};

template <int k> 
KdTree<k>::KdTree()
    : _sorted(0)
{}

template <int k>
KdTree<k>::~KdTree()
{}

template <int k>
void KdTree<k>::clear()
{
    _points.clear();
    _bbox.clear();
    _ids.clear();
    _sorted = 0;
}

template <int k>
void KdTree<k>::setPoints(float* p, int n)
{
    // copy points
    _points.resize(n);
    memcpy(&_points[0], p, sizeof(Point)*n);

    // compute bbox
    if (n) {
        _bbox.set(p);
        for (int i = 1; i < n; i++)
            _bbox.grow(_points[i].p);
    } else _bbox.clear();

    // assign sequential ids
    _ids.clear();
    _ids.reserve(n);
    while (int(_ids.size()) < n) _ids.push_back(_ids.size());
    _sorted = 0;
}

template <int k>
void KdTree<k>::allocatePoints(int n)
{
    _points.resize(n);
    _bbox.clear();
    
    _ids.clear();
    _sorted = 0;
}

template <int k>
void KdTree<k>::sort()
{
    if (_sorted) return;
    _sorted = 1;

    // reorder ids to sort points
    uint32_t np = _points.size();
    if (!np) return;
    while (_ids.size() < np) _ids.push_back(_ids.size());
    if (np > 1) sortSubtree(0, np, 0);

    // reorder points to match id order
    std::vector<Point> newpoints(np);
    for (uint32_t i = 0; i < np; i++)
        newpoints[i] = _points[_ids[i]];
    std::swap(_points, newpoints);
}

template <int k> 
void KdTree<k>::sortSubtree(int n, int size, int j)
{
    int left, right; ComputeSubtreeSizes(size, left, right);

    // partition range [n, n+size) along axis j into two subranges:
    //   [n, n+leftSize+1) and [n+leftSize+1, n+size)
    std::nth_element(&_ids[n], &_ids[n+left], &_ids[n+size], 
                     ComparePointsById(&_points[0].p[j]));
    // move median value (nth element) to front as root node of subtree
    std::swap(_ids[n], _ids[n+left]);

    // sort left and right subtrees using next discriminant 
    if (left <= 1) return;
    if (k > 1) j = (j+1)%k;
    sortSubtree(n+1, left, j);
    if (right <= 1) return;
    sortSubtree(n+left+1, right, j);
}

template <int k>
float KdTree<k>::findNPoints(std::vector<int>& result,
    std::vector<float>& distanceSquared,const float p[k],int nPoints,float maxRadius) const
{
    result.clear();
    distanceSquared.clear();
    float radius_squared=maxRadius*maxRadius;

    if (!size() || !_sorted || nPoints<1) return radius_squared;

    NearestNQuery query(result,distanceSquared,p,nPoints,radius_squared);
    findNPoints(query,0,size(),0);
    return query.maxRadiusSquared;
}

template<int k>
void KdTree<k>::findNPoints(KdTree<k>::NearestNQuery& query,int n,int size,int j) const
{
    const float* p = &_points[n].p[0];

    // Compute squared distance for this entry
    float pDistanceSquared=0;

    for (int axis=0;axis<k;axis++){
        float tmp=p[axis]-query.pquery[axis];
        pDistanceSquared+=tmp*tmp;
    }

    if (size>1) {
        // not at a leaf, visit children
        float axis_distance = query.pquery[j]-p[j];
        int left, right; ComputeSubtreeSizes(size, left, right);
        int nextj=(k>1)?(j+1)%k:j;

        if (axis_distance>0) { // visit right definitely, and left if within distance
            if (right) findNPoints(query, n+left+1, right,nextj);
            if ( axis_distance*axis_distance < query.maxRadiusSquared)
                findNPoints(query, n+1, left,nextj);
        } else { // visit left definitely, and right if within distance
            findNPoints(query, n+1, left,nextj);
            if (right && axis_distance*axis_distance < query.maxRadiusSquared)
                findNPoints(query, n+left+1, right,nextj);
        }
    }

    if (pDistanceSquared<query.maxRadiusSquared) {
        if (query.result.size()<query.maxPoints) {
            query.result.push_back(n);
            query.distanceSquared.push_back(pDistanceSquared);
            if (query.result.size()==query.maxPoints)
                query.maxRadiusSquared=buildHeap(query.result,query.distanceSquared);
        } else // already have heap, find somebody to throw out
            query.maxRadiusSquared=insertToHeap(query.result,query.distanceSquared,n,pDistanceSquared);
    }

}

template <int k>
int KdTree<k>::findNearest(const float p[k], float maxRadius) const
{
    if (!size() || !_sorted) return -1;

    NearestQuery query(p,maxRadius*maxRadius);
    findNearest(query,0,size(),0);
    return query.nearest;
}

template<int k>
void KdTree<k>::findNearest(KdTree<k>::NearestQuery& query,int n,int size,int j) const
{
    const float* p = &_points[n].p[0];

    // Compute squared distance for this entry
    float pDistanceSquared=0;

    for (int axis=0;axis<k;axis++){
        float tmp=p[axis]-query.pquery[axis];
        pDistanceSquared+=tmp*tmp;
    }

    if (pDistanceSquared<query.maxRadiusSquared) {
        query.nearest = n;
        query.maxRadiusSquared = pDistanceSquared;
    }

    if (size>1) {
        // not at a leaf, visit children
        float axis_distance = query.pquery[j]-p[j];
        int left, right; ComputeSubtreeSizes(size, left, right);
        int nextj=(k>1)?(j+1)%k:j;
        if (axis_distance>0) { // visit right definitely, and left if within distance
            if (right) findNearest(query, n+left+1, right,nextj);
            if ( axis_distance*axis_distance < query.maxRadiusSquared)
                findNearest(query, n+1, left,nextj);
        } else { // visit left definitely, and right if within distance
            findNearest(query, n+1, left,nextj);
            if (right && axis_distance*axis_distance < query.maxRadiusSquared)
                findNearest(query, n+left+1, right,nextj);
        }
    }
}

template <int k>
void KdTree<k>::findPoints(std::vector<int>& result, const BBox<k>& bbox) const
{
    if (!size() || !_sorted) return;
    if (!bbox.intersects(_bbox)) return;
    findPoints(result, bbox, 0, size(), 0);
}

template <int k>
void KdTree<k>::findPoints(std::vector<int>& result, const BBox<k>& bbox,
                           int n, int size, int j) const
{
    // check point at n for inclusion
    const float* p = &_points[n].p[0];
    if (bbox.inside(p))
        result.push_back(n);

    if (size == 1) return;

    // visit left subtree
    int left, right; ComputeSubtreeSizes(size, left, right);
    int nextj = (k > 1)? (j+1)%k : j;
    if (p[j] >= bbox.min[j])
        findPoints(result, bbox, n+1, left, nextj);

    // visit right subtree
    if (right && p[j] <= bbox.max[j])
        findPoints(result, bbox, n+left+1, right, nextj);
}

#endif
