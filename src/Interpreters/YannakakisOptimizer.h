#pragma once

#include <Interpreters/Aliases.h>
#include <Interpreters/DatabaseAndTableWithAlias.h>
#include <Interpreters/InDepthNodeVisitor.h>
#include <Storages/IStorage.h>
#include <Core/Settings.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>

namespace DB
{

class YannakakisOptimizer {

public:
    bool applyYannakakis(ASTSelectQuery & select);
};

void collectTablesAndPredicates(ASTSelectQuery & select, std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates);
void printTablesAndPredicates(const std::unordered_map<std::string, std::unordered_set<std::string>> tablesAndPredicates);

}

