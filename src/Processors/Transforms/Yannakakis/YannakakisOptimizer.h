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
    std::unordered_map<std::string, ASTPtr> &tableObjects,
    std::unordered_map<std::string, ASTPtr> &predicateObjects);
void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates);
void gyoReduction(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                  std::unordered_map<std::string, std::vector<std::string>> & joinTree,
                  DisjointSet &ds);
void addTableAndPredicateFromIdentifier(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                                        std::string identifier);
void bottomUpSemiJoin(ASTSelectQuery & select,
                      std::unordered_map<std::string, std::vector<std::string>> &join_tree,
                      std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                      std::unordered_map<std::string, ASTPtr> &tableObjects,
                      std::unordered_map<std::string, ASTPtr> &predicateObjects,
                      ASTPtr &subquery_template,
                      DisjointSet &ds);
void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> &joinTree, std::vector<std::string> &ordering);
void removeChild(ASTPtr & parent, const ASTPtr & child);


bool tryRemoveFromWhere(ASTSelectQuery & select, const ASTPtr & node);
bool bfsAnd(const ASTPtr & root, std::function<bool(ASTFunction &, const ASTPtr &)> visitor);
ASTPtr makeSubqueryTemplate(const String & table_alias);
ASTPtr buildJoinTreeRec(std::string &nodeIdentifier,
                        ASTSelectQuery &select,
                        std::unordered_map<std::string, ASTPtr> &tableObjects,
                        std::unordered_map<std::string, ASTPtr> &predicateObjects,
                        std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                        std::unordered_map<std::string, std::vector<std::string>> &adjacent,
                        ASTPtr  &subquery_template);
ASTs getJoinPredicates(std::string &leftIdentifier,
                       std::string &rightIdentifier,
                       std::unordered_map<std::string, ASTPtr> &predicateObjects,
                       std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                       DisjointSet &ds);
ASTPtr makeConjunction(const ASTs & nodes);
ASTPtr makeSubqueryQualifiedAsterisk();
void setTablesOfSubquery(ASTPtr & subquery, ASTs & tables);
bool isEquiJoin(ASTs functionArguments);
/*void removeChild(ASTPtr & parent, const ASTPtr & child);
*/
}

