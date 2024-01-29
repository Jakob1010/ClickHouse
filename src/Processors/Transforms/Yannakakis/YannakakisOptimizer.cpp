#include "YannakakisOptimizer.h"
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <Functions/FunctionsComparison.h>
#include <Functions/FunctionsLogical.h>
#include <Parsers/ASTFunction.h>
#include <Parsers/ASTQualifiedAsterisk.h>
#include <Parsers/ParserTablesInSelectQuery.h>
#include <Parsers/parseQuery.h>
#include "DisjointSet.h"
#include "Interpreters/IdentifierSemantic.h"
#include <Parsers/queryToString.h>
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"
#include <Parsers/ASTIdentifier.h>


namespace DB
{

namespace ErrorCodes
{
extern const int LOGICAL_ERROR;
}

/*
 * Assumptions that are made
 * - selection predicates in WHERE are always referred with a table alias e.g. table.attribute = 4
 *
 *
 */

bool YannakakisOptimizer::applyYannakakis(ASTSelectQuery & select)
{
    if (!select.tables())
        return false;
    if (select.tables()->as<ASTTablesInSelectQuery>()->children.size() < 2)
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;


    // stores tables and their predicates
    std::unordered_map<std::string, std::unordered_set<std::string>> tableWithPredicateNames;
    std::unordered_map<std::string, ASTPtr> tableObjects;
    std::unordered_map<std::string, ASTPtr> predicateObjects;
    std::unordered_map<std::string, std::vector<ASTPtr>> selectionObjects;
    DisjointSet ds; // every join predicate is represented by a EquivalenceClasses
                    // e.g.  all predicates (a.a, b.a, c.a) will be in the same class
                    // A INNER JOIN B ON(a.a == b.a) INNER JOIN C (b.a == c.a)

    // get tables and predicates
    collectTablesAndPredicates(select, tableWithPredicateNames, ds, tableObjects, predicateObjects, selectionObjects);
    printTablesAndPredicates(tableWithPredicateNames);

    // gyo reduction
    std::unordered_map<std::string, std::vector<std::string>> join_tree;
    if (gyoReduction(tableWithPredicateNames, join_tree, ds))
        removeJoin(select);

    // reroot based on selected attributes
    rerootTree(join_tree, "t");

    // apply Yannakakis algorithm
    bottomUpSemiJoin(select, join_tree, tableWithPredicateNames, tableObjects,
                     predicateObjects, ds, selectionObjects);

    ds.dumpSet();

    // Overwrite WHERE statement
    ASTPtr always_true = makeASTFunction("equals", std::make_shared<ASTLiteral>(1), std::make_shared<ASTLiteral>(1));
    select.setExpression(ASTSelectQuery::Expression::WHERE, std::move(always_true));

    std::cout << "Optimized Query" << std::endl;
    std::cout << select.dumpTree() << std::endl;
    std::cout << "applyYannakakis() finished" << std::endl;
    return true;
}

void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> & join_tree, std::vector<std::string> & ordering)
{
    if (join_tree.size() != 0)
        return;
    std::cout << "computeTopologicalOrdering()" << std::endl;
    std::unordered_map<std::string, int> incomingCnt; // holds count of incoming edges
    std::unordered_map<std::string, std::vector<std::string>> adjacent; // (node x) -> (directed edge) -> (node y)
    std::queue<std::string> q;

    // create incoming count and adjacent list
    for (const auto & kv : join_tree)
    {
        const std::string & edge = kv.first;
        const std::vector<std::string> & neighbors = kv.second;
        for (const std::string & neighbor : neighbors)
        {
            if (!incomingCnt.contains(edge))
                incomingCnt[edge] = 0;
            incomingCnt[edge] += 1;
            adjacent[neighbor].push_back(edge);
        }
    }

    // add nodes with no incoming edges to q
    for (const auto & kv : adjacent)
    {
        const std::string & node = kv.first;
        if (!incomingCnt.contains(node) || incomingCnt[node] == 0)
            q.push(node);
    }

    //  compute topological order
    while (!q.empty())
    {
        // Get and print the front element of the queue
        std::string child = q.front();
        q.pop();
        const std::vector<std::string> & neighbors = adjacent[child];
        for (const std::string & parent : neighbors)
        {
            ordering.push_back(child + "- SEMI JOIN -" + parent);
            incomingCnt[parent] -= 1;
            if (incomingCnt[parent] == 0)
                q.push(parent);
        }
    }
}

ASTPtr makeConjunction(const ASTs & nodes)
{
    if (nodes.empty())
        throw Exception(ErrorCodes::LOGICAL_ERROR, "Expected a non-zero number of nodes.");

    if (nodes.size() == 1)
        return nodes[0]->clone();

    ASTs arguments;
    arguments.reserve(nodes.size());
    for (const auto & ast : nodes)
        arguments.emplace_back(ast->clone());

    return makeASTFunction(NameAnd::name, std::move(arguments));
}

ASTs getJoinPredicates(
    std::string & leftTableIdentifier,
    std::string & rightTableIdentifier,
    std::string & rightJoinPredicateName,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    DisjointSet & ds)
{
    std::cout << "getJoinPredicates() for" << std::endl;
    std::cout << leftTableIdentifier << std::endl;
    std::cout << rightTableIdentifier << std::endl;

    ASTs joinPredicates;
    for (const auto & leftElement : tablesAndPredicates[leftTableIdentifier])
    {
        for (const auto & rightElement : tablesAndPredicates[rightTableIdentifier])
        {
            if (ds.find(leftTableIdentifier + "." + leftElement) == ds.find(rightTableIdentifier + "." + rightElement))
            {
                std::cout << "match" << std::endl;
                std::cout << leftElement << std::endl;
                std::cout << rightElement << std::endl;
                ASTPtr rightJoinPredicate = predicateObjects[rightTableIdentifier + "." + rightElement];
                ASTPtr leftJoinPredicate = predicateObjects[leftTableIdentifier + "." + leftElement];

                // Need to adapt join predicate name because it now comes from the subquery
                    if (rightTableIdentifier != rightJoinPredicateName) {
                        std::cout << "need to change tablename to subquery" << std::endl;
                        rightJoinPredicate = std::make_shared<ASTIdentifier>(std::vector<String>{rightJoinPredicateName, rightElement});
                    }
                auto equal_function = makeASTFunction(
                    "equals", leftJoinPredicate, rightJoinPredicate);
                joinPredicates.emplace_back(equal_function);
            }
        }
    }
    return joinPredicates;
}

ASTPtr buildJoinTreeRec(
    std::string & leftIdentifier,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, std::vector<std::string>> & adjacent,
    DisjointSet & ds,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates,
    int & subQueryCnt)
{
    std::cout << "buildJoinTreeRec()" << std::endl;
    std::cout << leftIdentifier << std::endl;
    const std::string subqueryName = "sub"+std::to_string(subQueryCnt);
    ASTPtr subquery = makeSubqueryTemplate(subqueryName);
    ASTs new_children;

    // add current root (left)
    ASTPtr leftTable = tableObjects[leftIdentifier];
    new_children.push_back(leftTable);

    ASTs whereFilter;
    addSelectionPredicatesOfTable(leftIdentifier, whereFilter, selectionPredicates);

    // add direct neighbors of "leftIdentifier"
    for (std::string & rightIdentifier : adjacent[leftIdentifier])
    {
        std::cout << rightIdentifier << std::endl;
        ASTPtr neighborPtr = tableObjects[rightIdentifier];
        ASTPtr right;
        std::string rightJoinPredicateName = rightIdentifier;
        std::cout << neighborPtr->as<ASTTablesInSelectQueryElement>()->dumpTree() << std::endl;

        // check if neighbor is leaf node
        if (!adjacent.contains(rightIdentifier) || adjacent[rightIdentifier].size() == 0)
        {
            right = tableObjects[rightIdentifier];
            addSelectionPredicatesOfTable(rightIdentifier, whereFilter, selectionPredicates);
        }
        else
        {
            subQueryCnt+=1;
            right = buildJoinTreeRec(
                rightIdentifier, tableObjects, predicateObjects, tablesAndPredicates,
                adjacent, ds, selectionPredicates, subQueryCnt);
            rightJoinPredicateName = "sub"+std::to_string(subQueryCnt);
        }
        auto * rightTable = right->as<ASTTablesInSelectQueryElement>();

        // create join
        auto join_ast = std::make_shared<ASTTableJoin>();
        ASTs join_predicates = getJoinPredicates(leftIdentifier, rightIdentifier, rightJoinPredicateName,
                                                 predicateObjects, tablesAndPredicates, ds);
        std::cout << "Number of Join Predicates" << std::endl;
        std::cout << join_predicates.size() << std::endl;
        if (join_predicates.empty())
            join_ast->kind = JoinKind::Cross;
        else {
            join_ast->kind = JoinKind::Left;
            join_ast->strictness = JoinStrictness::Semi;
            join_ast->on_expression = makeConjunction(join_predicates);
            join_ast->children.emplace_back(join_ast->on_expression);
        }
        // add join
        rightTable->table_join = join_ast;
        rightTable->children.emplace_back(join_ast);
        new_children.push_back(right);
    }

    ASTSelectQuery selectOfSubquery;
    setASTSelectQuery(subquery, selectOfSubquery);

    // set select
    auto expression_list = std::make_shared<ASTExpressionList>();
    expression_list->children.emplace_back(makeSubqueryQualifiedAsterisk(leftIdentifier));
    selectOfSubquery.setExpression(ASTSelectQuery::Expression::SELECT, std::move(expression_list));

    // set from
    setTablesOfSubquery(subquery, new_children);

    // set where
    if (whereFilter.size() >= 1)
        setSelectionOfSubquery(subquery, whereFilter);

    std::cout << "Subquery" << std::endl;
    std::cout << subquery->dumpTree() << std::endl;

    return subquery;
}

void addSelectionPredicatesOfTable(std::string & tableIdentifier, ASTs & whereSubquery,
                                   std::unordered_map<std::string, std::vector<ASTPtr>> & allPredicates) {
    if (allPredicates.contains(tableIdentifier))
        for (auto & wherePredicate : allPredicates[tableIdentifier])
            whereSubquery.push_back(wherePredicate);
}

void rerootTree(std::unordered_map<std::string, std::vector<std::string>> &join_tree, const std::string &newRoot) {
    std::cout << "rerootTree() with " << newRoot << std::endl;

    // join_tree only stores parent->children relationship,
    // but we also need children -> parent for reroot
    std::unordered_map<std::string, std::vector<std::string>> parents;

    // create parent relationship
    for (const auto & kv : join_tree)
    {
        const std::string & parent = kv.first;
        const std::vector<std::string> & children = kv.second;
        for (const std::string & child : children)
            parents[child].push_back(parent);
    }

    // Perform rerooting
    std::string current = newRoot;
    while (parents.contains(current)) { // newRoot is already root (has no parents)
        const std::string& p = parents[current][0]; // Each child has only one parent

        // remove current form child list of parent
        join_tree[p].erase(std::remove(join_tree[p].begin(), join_tree[p].end(), current), join_tree[p].end());

        // add parent as child from current
        join_tree[current].push_back(p);

        current = p;
    }

    std::cout << "Rerooted Join Tree" << std::endl;
    for (const auto & kv : join_tree)
    {
        const std::string & edge = kv.first;
        const std::vector<std::string> & neighbors = kv.second;
        std::cout << "Edge: " << edge << ", Neighbors: ";
        for (const std::string & neighbor : neighbors)
            std::cout << neighbor << " ";
        std::cout << std::endl;
    }

}

void setASTSelectQuery(ASTPtr & subquery, ASTSelectQuery & selectQuery) {
    if (!subquery)
        return;

    if (subquery->as<ASTSelectQuery>())
    {
        std::cout << "setASTSelectQuery()" << std::endl;
        selectQuery = *subquery->as<ASTSelectQuery>();
        return;
    }

    for (auto & childElement : subquery->children)
        setASTSelectQuery(childElement, selectQuery);

}

void setProjectionOfSubquery(ASTPtr & subquery) {

    if (!subquery)
        return;

    if (subquery->as<ASTSelectQuery>())
    {
        auto select_expr_list = std::make_shared<ASTExpressionList>();

        auto select_query = std::make_shared<ASTSelectQuery>();
        select_query->children.push_back(select_expr_list);

        select_query->setExpression(ASTSelectQuery::Expression::SELECT, select_expr_list);
        return;
    }

    for (auto & childElement : subquery->children)
        setProjectionOfSubquery(childElement);
}


void setSelectionOfSubquery(ASTPtr & subquery, ASTs & selectionPredicates)
{
    if (!subquery)
        return;

    if (subquery->as<ASTSelectQuery>())
    {
        std::cout << "set where" << std::endl;
        subquery->as<ASTSelectQuery>()->setExpression(ASTSelectQuery::Expression::WHERE, makeConjunction(selectionPredicates));
        return;
    }

    for (auto & childElement : subquery->children)
        setSelectionOfSubquery(childElement, selectionPredicates);
}

void setTablesOfSubquery(ASTPtr & subquery, ASTs & tables)
{
    if (!subquery)
        return;

    if (subquery->as<ASTTablesInSelectQuery>())
    {
        std::cout << "Found it bitch!" << std::endl;
        subquery->children = tables;
        return;
    }

    for (auto & childElement : subquery->children)
        setTablesOfSubquery(childElement, tables);
}


void bottomUpSemiJoin(
    ASTSelectQuery & select,
    std::unordered_map<std::string, std::vector<std::string>> & join_tree,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicateObjects,
    DisjointSet & ds,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates)
{
    auto * query_tables = select.tables()->as<ASTTablesInSelectQuery>();
    if (query_tables->children.size() < 2 || predicateObjects.size() < 1 || tableObjects.size() < 1)
        return;

    std::cout << "bottomUpSemiJoin()" << std::endl;
    std::unordered_map<std::string, int> incomingCnt; // holds count of incoming edges
    //std::unordered_map<std::string, std::vector<std::string>> adjacent; // (node x) -> (directed edge) -> (node y)
    std::string root;

    // create incoming count and adjacent list
    for (const auto & kv : join_tree)
    {
        const std::vector<std::string> & neighbors = kv.second;
        for (const std::string & neighbor : neighbors)
        {
            if (!incomingCnt.contains(neighbor))
                incomingCnt[neighbor] = 0;
            incomingCnt[neighbor] += 1;
            //adjacent[edge].push_back(neighbor);
        }
    }

    // find root
    for (const auto & kv : join_tree)
    {
        const std::string & node = kv.first;
        if (!incomingCnt.contains(node) || incomingCnt[node] == 0)
        {
            root = node;
            break;
        }
    }

    std::cout << "root" << std::endl;
    std::cout << root << std::endl;
    ASTPtr rootObject = tableObjects[root];
    int subQueryCnt = 0;
    ASTPtr subQuery = buildJoinTreeRec(root, tableObjects, predicateObjects, tablesAndPredicates,
                                       join_tree, ds, selectionPredicates, subQueryCnt);
    ASTs new_children;
    new_children.push_back(subQuery);
    select.tables()->as<ASTTablesInSelectQuery>()->children = new_children;

    // make asterix for subquery
    auto expression_list = std::make_shared<ASTExpressionList>();
    const std::string subqueryName = "sub"+std::to_string(0);
    expression_list->children.emplace_back(makeSubqueryQualifiedAsterisk(subqueryName));
    //select.setExpression(ASTSelectQuery::Expression::SELECT, std::move(expression_list));
}

ASTPtr makeSubqueryQualifiedAsterisk(const std::string & identifier)
{
    auto asterisk = std::make_shared<ASTQualifiedAsterisk>();
    asterisk->qualifier = std::make_shared<ASTIdentifier>(identifier);
    asterisk->children.push_back(asterisk->qualifier);
    return asterisk;
}

void addTableAndPredicateFromIdentifier(
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates, std::string identifier)
{
    size_t dotPos = identifier.find('.');
    if (dotPos != std::string::npos)
    {
        std::string table = identifier.substr(0, dotPos);
        std::string predicate = identifier.substr(dotPos + 1);
        tablesAndPredicates[table].insert(predicate);
    }
}

std::pair<std::string, std::string> splitIdentifier(const std::string & identifier)
{
    size_t dotPos = identifier.find('.');
    if (dotPos != std::string::npos)
    {
        std::string table = identifier.substr(0, dotPos);
        std::string predicate = identifier.substr(dotPos + 1);
        return {table, predicate};
    }
    // If there's no dot
    return {"", identifier};
}

bool isIdentifier(const ASTPtr & ast) {
    try
    {
        getIdentifierName(ast);
        return true;
    }
    catch (const DB::Exception & e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }
}

bool isEquiJoin(ASTs functionArguments)
{
    return isIdentifier(functionArguments[0]) && isIdentifier(functionArguments[1]);
}

std::string extractTableAliasAfterAS(const std::string& query) {
    std::size_t asPos = query.find(" AS ");
    if (asPos == std::string::npos) {
        return ""; // "AS" not found
    }

    std::size_t nameStart = asPos + 4; // 4 characters in " AS "
    std::size_t nameEnd = query.find(' ', nameStart);

    if (nameEnd == std::string::npos) {
        return query.substr(nameStart); // Extract till end if no space after name
    } else {
        return query.substr(nameStart, nameEnd - nameStart);
    }
}

void removeJoin(ASTSelectQuery & select) {
    std::cout << "removeJoins()" << std::endl;
    auto * ast_tables = select.tables()->as<ASTTablesInSelectQuery>();
    std::unordered_map<std::string, ASTPtr> tables;
    for (size_t i = 0; i < ast_tables->children.size(); ++i)
    {
        auto * table = ast_tables->children[i]->as<ASTTablesInSelectQueryElement>();
        if (!table or !table->table_join)
            continue;

        // remove old join
        if (table->table_join)
        {
            removeChild(ast_tables->children[i], table->table_join);
            table->table_join = nullptr;
        }
    }
}

void collectTablesAndPredicates(
    ASTSelectQuery & select,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    DisjointSet & ds,
    std::unordered_map<std::string, ASTPtr> & tableObjects,
    std::unordered_map<std::string, ASTPtr> & predicates,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates)
{
    std::cout << "collectTablesAndPredicates()" << std::endl;
    // join tables and predicates can be found in
    // FROM clause and WHERE clause

    // collect FROM
    std::cout << "collect FROM" << std::endl;
    auto * ast_tables = select.tables()->as<ASTTablesInSelectQuery>();
    std::unordered_map<std::string, ASTPtr> tables;
    for (size_t i = 0; i < ast_tables->children.size(); ++i)
    {
        auto * table = ast_tables->children[i]->as<ASTTablesInSelectQueryElement>();
        if (!table or !table->table_expression)
            continue;

        // set table object with table alias as key
        std::cout << "found table" << std::endl;
        std::string tableName = extractTableAliasAfterAS(queryToString(*(table)));
        std::cout << tableName << std::endl;
        tablesAndPredicates[tableName];
        tableObjects[tableName] = ast_tables->children[i];
    }

    // collect WHERE
    std::cout << "collect WHERE" << std::endl;
    if (select.where())
    {
        for (const std::shared_ptr<IAST> node : splitConjunctionsAst(select.where()))
        {
            if (const auto * func = node->as<ASTFunction>())
            {
                if (!func->arguments || func->arguments->children.size() != 2 || func->name != NameEquals::name)
                    continue;


                if (!isIdentifier(func->arguments->children[0]) || !isIdentifier(func->arguments->children[1])) {
                    // add selection predicate of type a=5
                    std::cout << node->dumpTree() << std::endl;
                    collectSelectionPredicate(node, func->arguments->children, tablesAndPredicates, selectionPredicates);
                    continue; // no join
                }

                const std::string identifierLeft = getIdentifierName(func->arguments->children[0]);
                const std::string identifierRight = getIdentifierName(func->arguments->children[1]);
                auto identifierLeftSplit = splitIdentifier(identifierLeft);
                auto identifierRightSplit = splitIdentifier(identifierRight);

                // make it an equivalent class if there is a join
                if (tablesAndPredicates.count(identifierLeftSplit.first) > 0 && tablesAndPredicates.count(identifierRightSplit.first) > 0)
                    ds.unionSets(identifierLeft, identifierRight);
                else {
                    // add selection predicate of type a="a"
                    std::cout << node->dumpTree() << std::endl;
                    collectSelectionPredicate(node, func->arguments->children, tablesAndPredicates, selectionPredicates);
                    continue; // no join
                }


                // store for every table identifier all of its join predicates
                // identifierSplit contains table as .first and predicate as .second
                tablesAndPredicates[identifierLeftSplit.first].insert(identifierLeftSplit.second);
                tablesAndPredicates[identifierRightSplit.first].insert(identifierRightSplit.second);

                // map table.predicate identifier to table.predicate
                predicates[identifierLeft] = func->arguments->children[0]->clone();
                predicates[identifierRight] = func->arguments->children[1]->clone();
            }
        }
    }
}

void collectSelectionPredicate(
    const ASTPtr & node,
    const ASTs & args,
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, std::vector<ASTPtr>> & selectionPredicates
    )
{
    std::cout << "collectSelectionPredicate()" << std::endl;
    if (selectionPredicates.size() > 100 || tablesAndPredicates.size() < 1 || args.size() < 1)
        return;

    bool addedPredicate = false;

    if(isIdentifier(args[0])) {
        const std::string identifierLeft = getIdentifierName(args[0]);
        auto identifierLeftSplit = splitIdentifier(identifierLeft);
        if (tablesAndPredicates.count(identifierLeftSplit.first) > 0) {
            addedPredicate = true;
            selectionPredicates[identifierLeftSplit.first].push_back(node);
        }
    }

    if(isIdentifier(args[1])) {
        const std::string identifierRight = getIdentifierName(args[1]);
        auto identifierRightSplit = splitIdentifier(identifierRight);
        if (tablesAndPredicates.count(identifierRightSplit.first) > 0) {
            addedPredicate = true;
            selectionPredicates[identifierRightSplit.first].push_back(node);
        }
    }

    if(!addedPredicate)
        std::cout << "collectSelectionPredicate -> predicate not added" << std::endl;
    else
        std::cout << "collectSelectionPredicate -> predicate added" << std::endl;
}


bool gyoReduction(
    std::unordered_map<std::string, std::unordered_set<std::string>> & tablesAndPredicates,
    std::unordered_map<std::string, std::vector<std::string>> & join_tree,
    DisjointSet & ds)
{
    std::cout << "gyoReduction()" << std::endl;

    // 3 steps:
    // i) deleting a vertex with degree 1 (i.e., a vertex occurring in a single edge)
    // ii) deleting an empty edge, or
    // iii) deleting an edge that is a subset of another edge

    // Create data structures to store relationships
    std::unordered_map<std::string, std::set<std::string>> vertex_to_edge;
    std::unordered_map<std::string, std::set<std::string>> edge_to_vertex;
    std::set<std::string> is_processed;
    bool success = false;

    // Read input data
    for (const auto & kv : tablesAndPredicates)
    {
        // edge of hypergraph is the table_name
        const std::string & edge = kv.first;
        // vertices of hypergraph are the attributes
        const std::unordered_set<std::string> & vertices = kv.second;

        // Associate vertices with hyperedges
        for (const std::string & vertex : vertices)
        {
            std::cout << vertex << std::endl;
            // equivalenceClass is always table_alias.attribute
            std::string equivalenceClass = edge + '.' + vertex;
            if (ds.exists(equivalenceClass)) {
                equivalenceClass = ds.find(equivalenceClass);
                std::cout << equivalenceClass << std::endl;
            }

            vertex_to_edge[equivalenceClass].insert(edge);
            edge_to_vertex[edge].insert(equivalenceClass);
        }
    }

    bool progress = true;
    while (edge_to_vertex.size() > 1 && progress)
    {
        progress = false;

        // i. Deleting vertices with degree 1
        auto cp_edge_to_vertex = edge_to_vertex;
        for (auto & kv : cp_edge_to_vertex)
        {
            const std::string & edge = kv.first;
            std::set<std::string> & vertices = kv.second;
            for (const std::string & vertex : vertices)
            {
                if (vertex_to_edge[vertex].size() == 1)
                { // Vertex degree == 1
                    edge_to_vertex[edge].erase(vertex);
                    vertex_to_edge[vertex].erase(edge);
                }
            }
        }

        // ii. & iii. Deleting empty edges or edges that are subsets of others
        cp_edge_to_vertex = edge_to_vertex;
        for (auto & kv : cp_edge_to_vertex)
        {
            const std::string & edge = kv.first;
            if (is_processed.count(edge) > 0)
                continue;
            bool is_a_subset = false;
            std::vector<std::string> subsets;
            for (auto & other_edge_kv : edge_to_vertex)
            {
                const std::string & other_edge = other_edge_kv.first;
                if (edge == other_edge)
                    continue;
                const std::set<std::string> & other_vertices = other_edge_kv.second;
                if (kv.second != other_vertices
                    && std::includes(other_vertices.begin(), other_vertices.end(), kv.second.begin(), kv.second.end()))
                {
                    is_a_subset = true;
                }
                else if (std::includes(kv.second.begin(), kv.second.end(), other_vertices.begin(), other_vertices.end()))
                {
                    subsets.push_back(other_edge);
                }
            }
            if (!is_a_subset)
            {
                progress = true;
                for (const std::string & other_edge : subsets)
                {
                    if (is_processed.count(other_edge) == 0)
                    {
                        is_processed.insert(other_edge);
                        join_tree[edge].push_back(other_edge);
                        // Update adjacency list
                        for (const std::string & vertex : edge_to_vertex[other_edge])
                            vertex_to_edge[vertex].erase(other_edge);
                        edge_to_vertex.erase(other_edge);
                    }
                }
            }
        }
    }

    if (edge_to_vertex.size() > 1) {
        std::cout << edge_to_vertex.size() << std::endl;
    }
    else
        success = true;

    // Print final state
    std::cout << "Join Tree" << std::endl;
    for (const auto & kv : join_tree)
    {
        const std::string & edge = kv.first;
        const std::vector<std::string> & neighbors = kv.second;
        std::cout << "Edge: " << edge << ", Neighbors: ";
        for (const std::string & neighbor : neighbors)
            std::cout << neighbor << " ";
        std::cout << std::endl;
    }

    std::cout << "gyo success" << std::endl;
    std::cout << success << std::endl;
    return success;
}

ASTPtr makeSubqueryTemplate(const String & table_alias)
{
    ParserTablesInSelectQueryElement parser(true);
    String query_template = "(select * from _t WHERE 1=1)";
    if (!table_alias.empty())
        query_template += " as " + table_alias;
    std::cout << "subquery was created" << std::endl;
    return parseQuery(parser, query_template, 0, DBMS_DEFAULT_MAX_PARSER_DEPTH);
}


void printTablesAndPredicates(std::unordered_map<std::string, std::unordered_set<std::string>> & tablePredicates)
{
    // Iterate through the tables and associated predicates
    for (const auto & pair : tablePredicates)
    {
        std::cout << "Table with Predicates" << std::endl;
        std::cout << "-" << pair.first << std::endl;
        for (const std::string & predicate : pair.second)
            std::cout << " --" << predicate << std::endl;
    }
}

void removeChild(ASTPtr & parent, const ASTPtr & child)
{
    if (child)
    {
        const auto * child_it
            = std::find_if(parent->children.begin(), parent->children.end(), [&](const auto & p) { return p.get() == child.get(); });
        if (child_it != parent->children.end())
            parent->children.erase(child_it);
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Couldn't find child to remove.");
    }
}


bool tryRemoveFromWhere(ASTSelectQuery & select, const ASTPtr & node)
{
    if (!node || !node->tryGetAlias().empty())
        return false;

    ASTFunction * rewritten_and = nullptr;
    /// First pass, convert AND(node, a) to AND(a)
    bfsAnd(
        select.refWhere(),
        [&](ASTFunction & parent, const ASTPtr & child)
        {
            if (child.get() == node.get())
            {
                if (!parent.tryGetAlias().empty())
                    return true;
                removeChild(parent.arguments, node);
                rewritten_and = &parent;
                return true;
            }
            return false; //continue traversal
        });

    if (!rewritten_and)
        return false;

    /// Second pass, convert AND(a) to a
    bool found = bfsAnd(
        select.refWhere(),
        [&](ASTFunction & parent, const ASTPtr & child)
        {
            if (child.get() != rewritten_and)
                return false; //continue traversal
            if (const auto * child_and = child->as<ASTFunction>();
                child_and && child_and->name == "and" && child_and->arguments->children.size() == 1)
            {
                auto * child_it
                    = std::find_if(parent.children.begin(), parent.children.end(), [&](const auto & p) { return p.get() == child_and; });
                if (child_it != parent.children.end())
                    *child_it = child_and->arguments->children.at(0);
            }
            return true;
        });
    if (!found)
    {
        /// If the AND(a) is not a child of the WHERE clause, it should *be* the WHERE clause. Rewrite to `WHERE a`.
        if (const auto * child = select.refWhere()->as<ASTFunction>(); child && child->name == "and")
        {
            if (child->arguments->children.size() == 1)
                select.refWhere() = child->arguments->children.at(0);
        }
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Rewrote an AND clause, but then lost track of it. This is a bug.");
    }
    return true;
}

bool bfsAnd(const ASTPtr & root, std::function<bool(ASTFunction &, const ASTPtr &)> visitor)
{
    ASTs parents = {root};

    for (size_t idx = 0; idx < parents.size();)
    {
        ASTPtr cur_parent = parents.at(idx);

        if (auto * function = cur_parent->as<ASTFunction>(); function && function->name == "and")
        {
            parents.erase(parents.begin() + idx);

            for (auto & child : function->arguments->children)
            {
                if (visitor(*function, child))
                    return true;
                parents.emplace_back(child);
            }
            continue;
        }
        ++idx;
    }
    return false;
}

}
