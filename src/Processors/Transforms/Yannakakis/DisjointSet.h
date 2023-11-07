#pragma once

#include <vector>
#include <unordered_map>
#include <string>

namespace DB
{
class DisjointSet
{
public:
    // Find the representative (root) of the set containing element x
    std::string find(std::string x);

    // Union two sets containing elements x and y
    void unionSets(std::string x, std::string y);

    // Prints key value pairs of parent
    void dumpSet();

    // Adds element to DS and returns true if successful
    bool addToSet(std::string);

    // Check if element exists
    bool exists(std::string);

private:
    std::unordered_map<std::string, std::string> parent;
    std::unordered_map<std::string, int> rank;
};

}
