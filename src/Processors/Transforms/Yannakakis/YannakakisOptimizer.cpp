#include "YannakakisOptimizer.h"
#include "DisjointSet.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "Functions/FunctionsLogical.h"
#include "Interpreters/IdentifierSemantic.h"
#include "Parsers/ASTIdentifier.h"
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"

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
    if(select.tables()->size()<=1)
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;

    // get tables and predicates
    std::unordered_map<std::string, std::unordered_set<std::string>> tablesWithPredicates;
    collectTablesAndPredicates(select, tablesWithPredicates);
    printTablesAndPredicates(tablesWithPredicates);

    // gyo reduction
    gyoReduction(tablesWithPredicates);

    std::cout << "applyYannakakis() finished" << std::endl;
    return true;
}

void collectTablesAndPredicates(ASTSelectQuery & select,
                                std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates)
{
    std::cout << "collectTablesAndPredicates()" << std::endl;
    const auto * ast_tables = select.tables()->as<ASTTablesInSelectQuery>();
    
    // collect predicates from where clause
    if (select.where()) {
        for (const auto& node : splitConjunctionsAst(select.where())) {
            // Check if the 'node' is not null
            if (node) {
                for (const auto& expressionList : node->children) {
                    // Check if the 'expressionList' is not null
                    if (expressionList) {
                        std::cout << expressionList->size() << std::endl;
                        for (const auto& identifier : expressionList->children) {
                            // Check if the 'identifier' is not null
                            // Assumption: we always use table_name.attribute
                            if (identifier) {
                                const std::string identifierName = getIdentifierName(identifier);
                                size_t dotPos = identifierName.find('.');
                                std::cout << identifierName << std::endl;
                                if (dotPos != std::string::npos) {
                                    std::string table = identifierName.substr(0, dotPos);
                                    std::string predicate = identifierName.substr(dotPos+1);
                                    tablesAndPredicates[table].insert(predicate);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // collect tables
    for (size_t i = 0; i < ast_tables->children.size(); ++i)
    {

        auto * table = ast_tables->children[i]->as<ASTTablesInSelectQueryElement>();
        if (!table or !table->table_expression)
            continue;
        for (const auto& identifier : table->table_expression->children) {
            const std::string tableName = getIdentifierName(identifier);
            tablesAndPredicates[tableName];
        }


        // extract predicates from join
        if (!table->table_join || !table->table_join->as<ASTTableJoin>())
            continue;



        /*
        auto * table_join = table->table_join->as<ASTTableJoin>();
        if (table_join)
        {
            std::cout << "Try to change join kind" << std::endl;
            table_join->kind = JoinKind::Left;
            table_join->strictness = JoinStrictness::Semi;

            if(table_join->on_expression){
             for (const auto & node : splitConjunctionsAst(table_join->on_expression))
            {
                    for(const auto & expressionList: node->children) {
                        for(const auto & identifier: expressionList->children) {
                            predicates.push_back(getIdentifierName(identifier));
                    }
                }
            }
             }

        }
         */
    }
}

void gyoReduction(std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates) {
    std::cout << "gyoReduction" << std::endl;
    // i) deleting a vertex with degree 1 (i.e., a vertex occurring in a single edge)
    // ii) deleting an empty edge, or
    // iii) deleting an edge that is a subset of another edge

    // Create data structures to store relationships
    std::unordered_map<std::string, std::set<std::string>> vertex_to_edge;
    std::unordered_map<std::string, std::set<std::string>> edge_to_vertex;
    std::map<std::string, std::vector<std::string>> join_tree;
    std::set<std::string> is_processed;

    // Read input data
    for (const auto& kv : tablesAndPredicates) {
        const std::string& edge = kv.first;
        const std::unordered_set<std::string>& vertices = kv.second;

        // Associate vertices with hyperedges
        for (const std::string& vertex : vertices) {
            vertex_to_edge[vertex].insert(edge);
            edge_to_vertex[edge].insert(vertex);
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

}
