#include <Functions/FunctionsLogical.h>
#include <Interpreters/YannakakisOptimizer.h>
#include <Interpreters/IdentifierSemantic.h>
#include <Parsers/ASTIdentifier.h>
#include <Parsers/ASTSelectQuery.h>
#include <Parsers/ASTTablesInSelectQuery.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>

namespace DB
{

namespace ErrorCodes
{
extern const int LOGICAL_ERROR;
}

using TablePredicatesMap = std::unordered_map<std::string, std::unordered_set<std::string>>;

bool YannakakisOptimizer::applyYannakakis(ASTSelectQuery & select)
{
    if(!select.tables())
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;

    // get tables and predicates
    std::unordered_map<std::string, std::unordered_set<std::string>> tablesWithPredicates;
    collectTablesAndPredicates(select, tablesWithPredicates);
    printTablesAndPredicates(tablesWithPredicates);

    std::cout << "applyYannakakis() finished" << std::endl;
    return true;
}

void collectTablesAndPredicates(ASTSelectQuery & select,
                                std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates)
{
    std::cout << "collectTablesAndPredicates()" << std::endl;
    const auto * ast_tables = select.tables()->as<ASTTablesInSelectQuery>();
    std::cout << ast_tables->children.size() << std::endl;

    // collect predicates from where clause
    if (select.where()) {
        for (const auto& node : splitConjunctionsAst(select.where())) {
            // Check if the 'node' is not null
            if (node) {
                for (const auto& expressionList : node->children) {
                    // Check if the 'expressionList' is not null
                    if (expressionList) {
                        for (const auto& identifier : expressionList->children) {
                            // Check if the 'identifier' is not null
                            // Assumption: we always use table_name.attribute
                            if (identifier) {
                                const std::string identifierName = getIdentifierName(identifier);
                                size_t dotPos = identifierName.find('.');
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

        /*
        // extract predicates from join
        if (!table->table_join || !table->table_join->as<ASTTableJoin>())
            continue;

        auto * table_join = table->table_join->as<ASTTableJoin>();
        if (table_join and table_join->on_expression)
        {
            for (const auto & node : splitConjunctionsAst(table_join->on_expression))
            {
                    for(const auto & expressionList: node->children) {
                        for(const auto & identifier: expressionList->children) {
                            predicates.push_back(getIdentifierName(identifier));
                    }
                }
            }
        }
        */
    }
}

void printTablesAndPredicates(const std::unordered_map<std::string, std::unordered_set<std::string>> tablePredicates) {
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
