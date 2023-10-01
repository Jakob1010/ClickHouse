#pragma once

#include <Interpreters/Aliases.h>
#include <Interpreters/DatabaseAndTableWithAlias.h>
#include <Interpreters/InDepthNodeVisitor.h>
#include <Storages/IStorage.h>
#include <Core/Settings.h>

namespace DB
{

class YannakakisOptimizer {

public:
    /*
    struct Data
    {
        const std::vector<StoragePtr> & storages;
        TablesWithColumns & tables; /// This will be reordered after visit() is called.
        const Aliases & aliases;
        const Settings & settings;
        const String current_database;
    };
     */

    //private:
    bool applyYannakakis(ASTSelectQuery & select, ASTPtr & ast, TablesWithColumns & tables);
};

}

