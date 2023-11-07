#include "YannakakisOptimizer.h"
#include "DisjointSet.h"
#include <queue>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include "Interpreters/IdentifierSemantic.h"
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"
#include <Parsers/ASTFunction.h>
#include <Functions/FunctionsComparison.h>
#include <Functions/FunctionsLogical.h>


namespace DB
{

namespace ErrorCodes
{
extern const int LOGICAL_ERROR;
}

bool YannakakisOptimizer::applyYannakakis(ASTSelectQuery & select)
{
    if(!select.tables())
        return false;
    if(select.tables()->size()<2)
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;

    // stores tables and their predicates
    std::unordered_map<std::string, std::unordered_set<std::string>> tableWithPredicateNames;
    std::unordered_map<std::string, ASTTablesInSelectQueryElement*> tableObjects;
    std::unordered_map<std::string, ASTPtr> predicateObjects;
    DisjointSet ds;     // every join predicate is represented by a EquivalenceClasses
                        // e.g.  all predicates (a.a, b.a, c.a) will be in the same class
                        // A INNER JOIN B ON(a.a == b.a) INNER JOIN C (b.a == c.a)

    // get tables and predicates
    collectTablesAndPredicates(select, tableWithPredicateNames, ds, tableObjects, predicateObjects);
    printTablesAndPredicates(tableWithPredicateNames);

    // gyo reduction
    std::unordered_map<std::string, std::vector<std::string>> join_tree;
    gyoReduction(tableWithPredicateNames, join_tree, ds);

    // get topological ordering
    std::vector<std::string> ordering;
    computeTopologicalOrdering(join_tree, ordering);

    // apply Yannakakis algorithm
    bottomUpSemiJoin(select, ordering);

    ds.dumpSet();
    std::cout << "Optimized Query" << std::endl;
    std::cout << select.dumpTree() << std::endl;
    std::cout << "applyYannakakis() finished" << std::endl;
    return true;
}

void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> & join_tree,
                                std::vector<std::string> ordering)
{
    if(join_tree.size() == 0)
        return;
    std::cout << "computeTopologicalOrdering()" << std::endl;
    std::unordered_map<std::string, int> incomingCnt; // holds count of incoming edges
    std::unordered_map<std::string, std::vector<std::string>> adjacent; // (node x) -> (directed edge) -> (node y)
    std::queue<std::string> q;

    // create incoming count and adjacent list
    for (const auto& kv : join_tree) {
        const std::string& edge = kv.first;
        const std::vector<std::string>& neighbors = kv.second;
        for (const std::string& neighbor : neighbors) {
            if(!incomingCnt.contains(edge))
                incomingCnt[edge] = 0;
            incomingCnt[edge] += 1;
            adjacent[neighbor].push_back(edge);
        }
    }

    // add nodes with no incoming edges to q
    for (const auto& kv : adjacent) {
        const std::string& node = kv.first;
        if(!incomingCnt.contains(node) || incomingCnt[node] == 0)
            q.push(node);
    }

    //  compute topological order
    while (!q.empty()) {
        // Get and print the front element of the queue
        std::string node = q.front();
        q.pop();
        const std::vector<std::string>& neighbors = adjacent[node];
        for (const std::string& neighbor: neighbors) {
            ordering.push_back(node + "->" + neighbor);
            incomingCnt[neighbor] -= 1;
            if (incomingCnt[neighbor] == 0) {
                q.push(neighbor);
            }
        }
    }

    std::cout << "Ordering of Joins" << std::endl;
    for (const std::string& node: ordering) {
        std::cout << node << std::endl;
    }

}

void bottomUpSemiJoin(ASTSelectQuery & select, std::vector<std::string> ordering)
{
    std::cout << "bottomUpSemiJoin()" << std::endl;
    auto * query_tables = select.tables()->as<ASTTablesInSelectQuery>();
    if (query_tables->children.size() != 100)
        return;

    /// Rules for query_tables->children (see ASTTablesInSelectQueryElement)
    /// First element either has el.table_expression or el.array_join
    /// The next elements either have both (table_join && table_expression) or (array_join)
    //  ASTs new_children;
    for (std::string & join : ordering)
    {
        std::cout << join << std::endl;
    }

}

void addTableAndPredicateFromIdentifier(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                                        std::string identifier)
{
        size_t dotPos = identifier.find('.');
        if (dotPos != std::string::npos)
        {
            std::string table = identifier.substr(0, dotPos);
            std::string predicate = identifier.substr(dotPos + 1);
            tablesAndPredicates[table].insert(predicate);
        }
}



void collectTablesAndPredicates(ASTSelectQuery & select,
                                std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                                DisjointSet & ds,
                                std::unordered_map<std::string, ASTTablesInSelectQueryElement*> &tableObjects,
                                std::unordered_map<std::string, ASTPtr> &predicates)
{
    std::cout << "collectTablesAndPredicates()" << std::endl;
    auto * ast_tables = select.tables()->as<ASTTablesInSelectQuery>();

    // collect predicates from where clause
    std::cout << "collectPredicates()" << std::endl;
    ASTs conjunctions;
    if (select.where()) {
        for (const std::shared_ptr<IAST> node: splitConjunctionsAst(select.where())) {
                if (const auto * func = node->as<ASTFunction>())
                {
                    if (!func->arguments || func->arguments->children.size() != 2 || func->name != NameEquals::name)
                        continue;
                    const std::string identifierName1 = getIdentifierName(func->arguments->children[0]);
                    const std::string identifierName2 = getIdentifierName(func->arguments->children[1]);
                    addTableAndPredicateFromIdentifier(tablesAndPredicates, identifierName1);
                    addTableAndPredicateFromIdentifier(tablesAndPredicates, identifierName2);
                    ds.unionSets(identifierName1, identifierName2);
                    predicates[ds.find(identifierName1)] = func->arguments->children[0]->clone();
                    predicates[ds.find(identifierName2)] = func->arguments->children[1]->clone();
                    conjunctions.push_back(node);
                }
        }
    }

    // collect tables
    std::cout << "collectTables()" << std::endl;
    for (size_t i = 0; i < ast_tables->children.size(); ++i)
    {
        auto * table = ast_tables->children[i]->as<ASTTablesInSelectQueryElement>();
        if (!table or !table->table_expression)
            continue;

        for (const auto& identifier : table->table_expression->children) {
            const std::string tableName = getIdentifierName(identifier);
            tablesAndPredicates[tableName];
            tableObjects[tableName] = table;
        }

        if (table->table_join && table->table_join->as<ASTTableJoin>())
        {
            std::cout << "Change JoinStrictness" << std::endl;
            auto join_ast = std::make_shared<ASTTableJoin>();
            join_ast->kind = JoinKind::Left;
            join_ast->strictness = JoinStrictness::Semi;
            if (i <= conjunctions.size()) {
                    std::cout << "add on expression" << std::endl;
                    join_ast->on_expression = conjunctions[i-1]->clone();
                    join_ast->children.emplace_back(join_ast->on_expression);
                    removeChild(ast_tables->children[i], table->table_join);
                    table->table_join = join_ast;
                    table->children.emplace_back(join_ast);
            }

        }

    }
}

void gyoReduction(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
                  std::unordered_map<std::string, std::vector<std::string>> & join_tree,
                  DisjointSet &ds) {
    std::cout << "gyoReduction()" << std::endl;
    // 3 steps:
    // i) deleting a vertex with degree 1 (i.e., a vertex occurring in a single edge)
    // ii) deleting an empty edge, or
    // iii) deleting an edge that is a subset of another edge

    // Create data structures to store relationships
    std::unordered_map<std::string, std::set<std::string>> vertex_to_edge;
    std::unordered_map<std::string, std::set<std::string>> edge_to_vertex;
    std::set<std::string> is_processed;

    // Read input data
    for (const auto& kv : tablesAndPredicates) {
        const std::string& edge = kv.first;
        const std::unordered_set<std::string>& vertices = kv.second;

        // Associate vertices with hyperedges
        for (const std::string& vertex : vertices) {
            std::string equivalenceClass = vertex;
            if(ds.exists(equivalenceClass))
                    equivalenceClass = ds.find(equivalenceClass);
            vertex_to_edge[equivalenceClass].insert(edge);
            edge_to_vertex[edge].insert(equivalenceClass);
        }
    }

    bool progress = true;
    while (edge_to_vertex.size() > 1 && progress) {
        progress = false;

        // i. Deleting vertices with degree 1
        auto cp = edge_to_vertex;
        for (auto& kv : cp) {
            const std::string& edge = kv.first;
            std::set<std::string>& vertices = kv.second;
            for (const std::string& vertex : vertices) {
                if (vertex_to_edge[vertex].size() == 1) { // Vertex degree == 1
                    edge_to_vertex[edge].erase(vertex);
                    vertex_to_edge[vertex].erase(edge);
                }
            }
        }

        // ii. & iii. Deleting empty edges or edges that are subsets of others
        cp = edge_to_vertex;
        for (auto& kv : cp) {
            const std::string& edge = kv.first;
            if (is_processed.count(edge) > 0) {
                continue;
            }
            bool is_a_subset = false;
            std::vector<std::string> subsets;
            for (auto& other_edge_kv : edge_to_vertex) {
                const std::string& other_edge = other_edge_kv.first;
                if (edge == other_edge) {
                    continue;
                }
                const std::set<std::string>& other_vertices = other_edge_kv.second;
                if (kv.second != other_vertices &&
                    std::includes(other_vertices.begin(), other_vertices.end(), kv.second.begin(), kv.second.end())) {
                    is_a_subset = true;
                } else {
                    subsets.push_back(other_edge);
                }
            }
            if (!is_a_subset) {
                progress = true;
                for (const std::string& other_edge : subsets) {
                    if (is_processed.count(other_edge) == 0) {
                        is_processed.insert(other_edge);
                        join_tree[edge].push_back(other_edge);
                        // Update adjacency list
                        for (const std::string& vertex : edge_to_vertex[other_edge]) {
                            vertex_to_edge[vertex].erase(other_edge);
                        }
                        edge_to_vertex.erase(other_edge);
                    }
                }
            }
        }
    }

    // Print final state
    std::cout << "Join Tree" << std::endl;
    for (const auto& kv : join_tree) {
        const std::string& edge = kv.first;
        const std::vector<std::string>& neighbors = kv.second;
        std::cout << "Edge: " << edge << ", Neighbors: ";
        for (const std::string& neighbor : neighbors) {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }
}


void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> &tablePredicates) {
    // Iterate through the tables and associated predicates
    for (const auto& pair : tablePredicates) {
        std::cout << "Table with Predicates"<< std::endl;
        std::cout << "-" << pair.first << std::endl;
        for (const std::string& predicate : pair.second) {
            std::cout << " --" << predicate << std::endl;
        }
    }
}

void removeChild(ASTPtr & parent, const ASTPtr & child)
{
    if (child)
    {
        const auto * child_it = std::find_if(parent->children.begin(), parent->children.end(), [&](const auto & p)
                                             {
                                                 return p.get() == child.get();
                                             });
        if (child_it != parent->children.end())
            parent->children.erase(child_it);
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Couldn't find child to remove.");
    }
}

}
