#include "YannakakisOptimizer.h"
#include "DisjointSet.h"
#include <queue>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <set>
#include "Interpreters/IdentifierSemantic.h"
#include "Parsers/ASTSelectQuery.h"
#include "Parsers/ASTTablesInSelectQuery.h"
#include "Parsers/ASTAsterisk.h"
#include <Parsers/parseQuery.h>
#include <Parsers/ASTFunction.h>
#include <Parsers/ASTQualifiedAsterisk.h>
#include <Parsers/ParserTablesInSelectQuery.h>
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
    if(select.tables()->as<ASTTablesInSelectQuery>()->children.size()<3)
        return false;

    std::cout << "applyYannakakis()" << std::endl;
    std::cout << select.dumpTree() << std::endl;

    // stores tables and their predicates
    std::unordered_map<std::string, std::unordered_set<std::string>> tableWithPredicateNames;
    std::unordered_map<std::string, ASTPtr> tableObjects;
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
    ParserTablesInSelectQueryElement parser(true);
    static ASTPtr subquery_template = parseQuery(parser, "(select * from _t) as `--.s`", 0, DBMS_DEFAULT_MAX_PARSER_DEPTH);
    std::cout << "Subquery Template" << std::endl;
    std::cout << subquery_template->dumpTree() << std::endl;
    bottomUpSemiJoin(select, join_tree, tableWithPredicateNames,
                     tableObjects, predicateObjects, subquery_template, ds);

    ds.dumpSet();

    // Overwrite WHERE statement
    ASTPtr always_true = makeASTFunction("equals", std::make_shared<ASTLiteral>(1), std::make_shared<ASTLiteral>(1));
    select.setExpression(ASTSelectQuery::Expression::WHERE, std::move(always_true));

    std::cout << "Optimized Query" << std::endl;
    std::cout << select.dumpTree() << std::endl;
    std::cout << "applyYannakakis() finished" << std::endl;
    return true;
}

void computeTopologicalOrdering(std::unordered_map<std::string, std::vector<std::string>> & join_tree,
                                std::vector<std::string> & ordering)
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
        std::string child = q.front();
        q.pop();
        const std::vector<std::string>& neighbors = adjacent[child];
        for (const std::string& parent: neighbors) {
            ordering.push_back(child + "- SEMI JOIN -" + parent);
            incomingCnt[parent] -= 1;
            if (incomingCnt[parent] == 0) {
                q.push(parent);
            }
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

ASTs getJoinPredicates(std::string &leftIdentifier,
                       std::string &rightIdentifier,
                       std::unordered_map<std::string, ASTPtr> &predicateObjects,
                       std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                       DisjointSet &ds)
{
    std::cout << "getJoinPredicates" << std::endl;
    ASTs joinPredicates;
    for (const auto &leftElement : tablesAndPredicates[leftIdentifier]) {
        for (const auto &rightElement: tablesAndPredicates[rightIdentifier]) {
            std::cout << "Compare" << std::endl;
            std::cout << leftElement << std::endl;
            std::cout << rightElement << std::endl;
            if (ds.find(leftElement) == ds.find(rightElement)) {
                auto equal_function =
                    makeASTFunction("equals", predicateObjects[leftIdentifier + "." + leftElement],
                                    predicateObjects[rightIdentifier + "." + rightElement]);
                joinPredicates.emplace_back(equal_function);
            }
        }
    }
    return joinPredicates;
}

ASTPtr buildJoinTreeRec(std::string &leftIdentifier,
                        ASTSelectQuery &select,
                        std::unordered_map<std::string, ASTPtr> &tableObjects,
                        std::unordered_map<std::string, ASTPtr> &predicateObjects,
                        std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                        std::unordered_map<std::string, std::vector<std::string>> &adjacent,
                        ASTPtr  &subquery_template,
                        DisjointSet &ds)
{
    std::cout << "buildJoinTreeRec()" << std::endl;
    std::cout << leftIdentifier << std::endl;
    ASTPtr subquery = subquery_template->clone();
    ASTPtr leftTable = tableObjects[leftIdentifier];
    ASTs new_children;
    new_children.push_back(leftTable);

    for(std::string &rightIdentifier : adjacent[leftIdentifier]) {
        std::cout << rightIdentifier << std::endl;
        ASTPtr neighborPtr = tableObjects[rightIdentifier];
        ASTPtr right;
        std::cout << neighborPtr->as<ASTTablesInSelectQueryElement>()->dumpTree() << std::endl;

        // check if neighbor is leaf node
        if(!adjacent.contains(rightIdentifier) || adjacent[rightIdentifier].size() == 0) {
            right = tableObjects[rightIdentifier];
        } else {
            right = buildJoinTreeRec(leftIdentifier, select, tableObjects, predicateObjects, tablesAndPredicates, adjacent, subquery_template, ds);
        }
        auto * rightTable = right->as<ASTTablesInSelectQueryElement>();

        // create join
        auto join_ast = std::make_shared<ASTTableJoin>();
        ASTs join_predicates = getJoinPredicates(leftIdentifier, rightIdentifier, predicateObjects, tablesAndPredicates, ds);
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

    setTablesOfSubquery(subquery, new_children);
    std::cout << "Subquery" << std::endl;
    std::cout << subquery->dumpTree() << std::endl;

    return subquery;
}

void setTablesOfSubquery(ASTPtr & subquery, ASTs & tables) {
    if (!subquery)
        return;

    if(subquery->as<ASTTablesInSelectQuery>()) {
        std::cout << "Found it bitch!" << std::endl;
        subquery->children = tables;
    }

    for (auto & childElement: subquery->children) {
        setTablesOfSubquery(childElement, tables);
    }
}



void bottomUpSemiJoin(ASTSelectQuery & select,
                      std::unordered_map<std::string, std::vector<std::string>> &join_tree,
                      std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                      std::unordered_map<std::string, ASTPtr> &tableObjects,
                      std::unordered_map<std::string, ASTPtr> &predicateObjects,
                      ASTPtr  &subquery_template,
                      DisjointSet &ds)
{
    auto * query_tables = select.tables()->as<ASTTablesInSelectQuery>();
    if (query_tables->children.size() < 2 || predicateObjects.size() < 1 || tableObjects.size() < 1)
        return;

    std::cout << "bottomUpSemiJoin()" << std::endl;
    std::unordered_map<std::string, int> incomingCnt; // holds count of incoming edges
    std::unordered_map<std::string, std::vector<std::string>> adjacent; // (node x) -> (directed edge) -> (node y)
    std::string root;

    // create incoming count and adjacent list
    for (const auto& kv : join_tree) {
        const std::string& edge = kv.first;
        const std::vector<std::string>& neighbors = kv.second;
        for (const std::string& neighbor : neighbors) {
            if(!incomingCnt.contains(neighbor))
                incomingCnt[neighbor] = 0;
            incomingCnt[neighbor] += 1;
            adjacent[edge].push_back(neighbor);
        }
    }

    // find root
    for (const auto& kv : adjacent) {
        const std::string& node = kv.first;
        if(!incomingCnt.contains(node) || incomingCnt[node] == 0) {
            root = node;
            break;
        }
    }

    std::cout << "root" << std::endl;
    std::cout << root << std::endl;
    ASTPtr rootObject = tableObjects[root];
    ASTPtr subQuery = buildJoinTreeRec(root, select, tableObjects, predicateObjects, tablesAndPredicates, adjacent, subquery_template, ds);
    ASTs new_children;
    new_children.push_back(subQuery);
    select.tables()->as<ASTTablesInSelectQuery>()->children = new_children;

    // make asterix for subquery
    auto expression_list = std::make_shared<ASTExpressionList>();
    expression_list->children.emplace_back(makeSubqueryQualifiedAsterisk());
    select.setExpression(ASTSelectQuery::Expression::SELECT, std::move(expression_list));
}

ASTPtr makeSubqueryQualifiedAsterisk()
{
    auto asterisk = std::make_shared<ASTQualifiedAsterisk>();
    asterisk->qualifier = std::make_shared<ASTIdentifier>("--.s");
    asterisk->children.push_back(asterisk->qualifier);
    return asterisk;
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

std::pair<std::string, std::string> splitIdentifier(const std::string& identifier)
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

bool isEquiJoin(ASTs functionArguments) {
    try {
        getIdentifierName(functionArguments[0]);
        getIdentifierName(functionArguments[1]);
        return true;
    } catch(const DB::Exception& e) {
        std::cout << e.what() << std::endl;
        return false;
    }
}

void collectTablesAndPredicates(ASTSelectQuery & select,
                                std::unordered_map<std::string, std::unordered_set<std::string>> &tablesAndPredicates,
                                DisjointSet & ds,
                                std::unordered_map<std::string, ASTPtr> &tableObjects,
                                std::unordered_map<std::string, ASTPtr> &predicates)
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

        if (table->table_join) {
            removeChild(ast_tables->children[i], table->table_join);
            table->table_join = nullptr;
        }

        for (const auto& identifier : table->table_expression->children) {
            const std::string tableName = getIdentifierName(identifier);
            tablesAndPredicates[tableName];
            tableObjects[tableName] = ast_tables->children[i];
        }
    }

    // collect WHERE
    std::cout << "collect WHERE" << std::endl;
    if (select.where()) {
        for (const std::shared_ptr<IAST> node: splitConjunctionsAst(select.where())) {
                if (const auto * func = node->as<ASTFunction>())
                {
                    if (!func->arguments || func->arguments->children.size() != 2 || func->name != NameEquals::name)
                        continue;

                    if (!isEquiJoin(func->arguments->children))
                        continue; // TODO: add other predicates like a=5
                    else
                        std::cout << "add predicate" << std::endl;

                    const std::string identifierLeft = getIdentifierName(func->arguments->children[0]);
                    const std::string identifierRight = getIdentifierName(func->arguments->children[1]);
                    auto identifierLeftSplit = splitIdentifier(identifierLeft);
                    auto identifierRightSplit = splitIdentifier(identifierRight);

                    // make it an equivalent class if there is a join
                    if (tablesAndPredicates.count(identifierLeftSplit.first) > 0
                        && tablesAndPredicates.count(identifierRightSplit.first) > 0)
                        ds.unionSets(identifierLeft, identifierRight);
                    else
                        continue; // TODO: add other predicates like a="a"

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

ASTPtr makeSubqueryTemplate(const String & table_alias)
{
    ParserTablesInSelectQueryElement parser(true);
    String query_template = "(select * from _t)";
    if (!table_alias.empty())
        query_template += " as " + table_alias;
    ASTPtr subquery_template = parseQuery(parser, "(select * from _t) as --.s", 100000000, 1000);
    if (!subquery_template)
        throw Exception(ErrorCodes::LOGICAL_ERROR, "Cannot parse subquery template");

    return subquery_template;

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


bool tryRemoveFromWhere(ASTSelectQuery & select, const ASTPtr & node)
{
    if (!node || !node->tryGetAlias().empty())
        return false;

    ASTFunction * rewritten_and = nullptr;
    /// First pass, convert AND(node, a) to AND(a)
    bfsAnd(select.refWhere(), [&](ASTFunction & parent, const ASTPtr & child)
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
    bool found = bfsAnd(select.refWhere(), [&](ASTFunction & parent, const ASTPtr & child)
                        {
                            if (child.get() != rewritten_and)
                                return false; //continue traversal
                            if (const auto * child_and = child->as<ASTFunction>(); child_and && child_and->name == "and" && child_and->arguments->children.size() == 1)
                            {
                                auto * child_it = std::find_if(parent.children.begin(), parent.children.end(), [&](const auto & p)
                                                               {
                                                                   return p.get() == child_and;
                                                               });
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
            {
                select.refWhere() = child->arguments->children.at(0);
            }
        }
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Rewrote an AND clause, but then lost track of it. This is a bug.");
    }
    return true;
}

bool bfsAnd(const ASTPtr & root, std::function<bool(ASTFunction &, const ASTPtr &)> visitor)
{
    ASTs parents = { root };

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

