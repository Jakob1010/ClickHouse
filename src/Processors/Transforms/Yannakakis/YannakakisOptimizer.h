#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "DisjointSet.h"
#include "Core/Settings.h"
#include "Interpreters/Aliases.h"
#include "Interpreters/DatabaseAndTableWithAlias.h"
#include "Interpreters/InDepthNodeVisitor.h"
#include "Storages/IStorage.h"

namespace DB
{

class YannakakisOptimizer {

public:
    bool applyYannakakis(ASTSelectQuery & select);
};

void collectTablesAndPredicates(
    ASTSelectQuery & select,
    std::unordered_map<std::string,std::unordered_set<std::string>> &tablesAndPredicates);
void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates);
void gyoReduction(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates);
}

