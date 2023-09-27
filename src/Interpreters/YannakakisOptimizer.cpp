#include <Functions/FunctionsLogical.h>
#include <Interpreters/YannakakisOptimizer.h>
#include <Interpreters/IdentifierSemantic.h>
#include <Parsers/ASTIdentifier.h>
#include <Parsers/ASTSelectQuery.h>
#include <Parsers/ASTTablesInSelectQuery.h>
#include <Storages/StorageJoin.h>
#include <iostream>

namespace DB
{

namespace ErrorCodes
{
    extern const int LOGICAL_ERROR;
}

bool YannakakisOptimizer::applyYannakakis(ASTSelectQuery & select, ASTPtr & query, TablesWithColumns & tables)
{
    if(!select.tables() or !query or tables.size() < 2)
        return false;

    const auto * query_tables = select.tables()->as<ASTTablesInSelectQuery>();
    if (!query_tables || query_tables->children.size() <= 1)
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;

    return true;
}
}
