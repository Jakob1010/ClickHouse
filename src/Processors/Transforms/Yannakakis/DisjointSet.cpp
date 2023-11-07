#include "DisjointSet.h"
#include <iostream>
namespace DB
{

bool DisjointSet::exists(std::string key)
{
    return parent.contains(key);
}

bool DisjointSet::addToSet(std::string key)
{
    if (parent.contains(key))
        return false;

    parent[key] = key;
    rank[key] = 1;

    return true;
}

// Find the representative (root) of the set containing element x
std::string DisjointSet::find(std::string key)
{
    if (!parent.contains(key))
        return "";

    if (key != parent[key])
        parent[key] = find(parent[key]);

    return parent[key];
}

// Union two sets containing elements x and y
void DisjointSet::unionSets(std::string x, std::string y)
{
    if (!parent.contains(x))
        addToSet(x);
    if(!parent.contains(y))
        addToSet(y);
    std::string rootX = find(x);
    std::string rootY = find(y);

    if (rootX != rootY)
    {
        // Union by rank to keep the tree balanced
        if (rank[rootX] < rank[rootY])
        {
            parent[rootX] = rootY;
        }
        else if (rank[rootX] > rank[rootY])
        {
            parent[rootY] = rootX;
        }
        else
        {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
    }
}

// Function to count and print parent key-value pairs
void DisjointSet::dumpSet()
{
    std::cout << "DisjointSet:" << std::endl;
    for (const auto& pair : parent)
        {
        std::cout << "Node: " << pair.first << ", Parent: " << pair.second << std::endl;
        }
}

}
