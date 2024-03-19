#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <Parsers/ASTFunction.h>
#include "Core/Settings.h"
#include "DisjointSet.h"
#include "Interpreters/Aliases.h"
#include "Interpreters/DatabaseAndTableWithAlias.h"
#include "Interpreters/IdentifierSemantic.h"
#include "Interpreters/InDepthNodeVisitor.h"
#include "Parsers/ASTIdentifier.h"
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"
#include "Storages/IStorage.h"

namespace DB
{

class YannakakisOptimizer
{
public:
    bool applyYannakakis(ASTSelectQuery & select);
};

void collectTablesAndPredicates(
    ASTSelectQuery & select,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    DisjointSet & ds,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicates,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates);
void collectSelectionPredicate(const ASTPtr & node,
                               const ASTs & args,
                               std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                               std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates
                               );
void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates);
bool gyoReduction(
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, std::vector<std::string>> & joinTree,
    DisjointSet & ds);
void addTableAndPredicateFromIdentifier(
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates, std::string identifier);
void bottomUpSemiJoin(
    ASTSelectQuery & select,
    std::unordered_map<std::string, std::vector<std::string>> & join_tree,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    DisjointSet & ds,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates);
void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> & joinTree, std::vector<std::string> & ordering);
void removeChild(ASTPtr & parent, const ASTPtr & child);


bool tryRemoveFromWhere(ASTSelectQuery & select, const ASTPtr & node);
bool bfsAnd(const ASTPtr & root, std::function<bool(ASTFunction &, const ASTPtr &)> visitor);
ASTPtr makeSubqueryTemplate(const String & table_alias);
ASTPtr buildJoinTreeRec(
    std::string & nodeIdentifier,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, std::vector<std::string>> & adjacent,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates,
    int & subQueryCnt);
ASTs getJoinPredicates(
    std::string & leftIdentifier,
    std::string & rightIdentifier,
    std::string & rightJoinPredicateName,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    DisjointSet & ds);
ASTPtr makeConjunction(const ASTs & nodes);
ASTPtr makeSubqueryQualifiedAsterisk(const std::string & identifier);
void addSelectionPredicatesOfTable(std::string & tableIdentifier, ASTs & whereSubquery,
                                   std::unordered_map<std::string, std::vector<ASTPtr>> & allPredicates);
void setTablesOfSubquery(ASTPtr & subquery, ASTs & tables);
void setSelectionOfSubquery(ASTPtr & subquery, ASTs & selectionPredicates);
void removeJoin(ASTSelectQuery & select);
bool isEquiJoin(ASTs functionArguments);
bool isIdentifier(const ASTPtr & ast);
std::string extractTableAliasAfterAS(const std::string& input);
bool rerootTree(std::unordered_map<std::string, std::vector<std::string>> &join_tree, ASTSelectQuery & selectQuery);
void setASTSelectQuery(ASTPtr & subquery, ASTSelectQuery & selectQuery);
/*void removeChild(ASTPtr & parent, const ASTPtr & child);
*/
}
