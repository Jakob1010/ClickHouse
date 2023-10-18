#pragma once

#include <vector>

class DisjointSet
{
public:
    DisjointSet(int size);

    // Find the representative (root) of the set containing element x
    int find(int x);

    // Union two sets containing elements x and y
    void unionSets(int x, int y);

    // Check if elements x and y are in the same set
    bool sameSet(int x, int y);

private:
    std::vector<int> parent;
    std::vector<int> rank;
};
