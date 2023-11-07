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
#include "Interpreters/IdentifierSemantic.h"
#include "Parsers/ASTIdentifier.h"
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"
#include <Parsers/ASTFunction.h>

namespace DB
{

class YannakakisOptimizer {

public:
    bool applyYannakakis(ASTSelectQuery & select);
};

void collectTablesAndPredicates(
    ASTSelectQuery & select,
    std::unordered_map<std::string,std::unordered_set<std::string>> &tablesAndPredicates,
    DisjointSet &ds,
    std::unordered_map<std::string, ASTTablesInSelectQueryElement*> &tableObjects,
    std::unordered_map<std::string, ASTPtr> &predicates);
void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates);
void gyoReduction(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                  std::unordered_map<std::string, std::vector<std::string>> & joinTree,
                  DisjointSet &ds);
void addTableAndPredicateFromIdentifier(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                                        std::string identifier);
void bottomUpSemiJoin(ASTSelectQuery & select, std::vector<std::string> ordering);
void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> &joinTree, std::vector<std::string> ordering);
void removeChild(ASTPtr & parent, const ASTPtr & child);
/*
bool tryRemoveFromWhere(ASTSelectQuery & select, const ASTPtr & node);
bool bfsAnd(const ASTPtr & root, std::function<bool(ASTFunction &, const ASTPtr &)> visitor);
void removeChild(ASTPtr & parent, const ASTPtr & child);
*/
}

