#include "interface/plugin.hpp"
#include "interface/database.hpp"
#include "interface/histogram.hpp"
#include "interface/variables.hpp"

#include "test/utils.hpp"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace theta;
using namespace theta::plugin;

BOOST_AUTO_TEST_SUITE(sqlite)


BOOST_AUTO_TEST_CASE(largefile){
   int argc = boost::unit_test::framework::master_test_suite().argc;
   char ** argv = boost::unit_test::framework::master_test_suite().argv;
   bool test = false;
   for(int i=1; i<argc; ++i){
       if(argv[i] == string("--sqlite_largefile")) test = true;
   }
   if(!test) return;
   BOOST_TEST_CHECKPOINT("enter sqlite largefile test");
   load_core_plugins();
   string config = "type = \"sqlite_database\"; filename = \"test_largefile.db\";";
   boost::shared_ptr<VarIdManager> vm(new VarIdManager());
   ConfigCreator cc(config, vm);
   Configuration cfg = cc.get();
   std::auto_ptr<Database> db = PluginManager<Database>::instance().build(cfg);
   std::auto_ptr<Table> table = db->create_table("test_table");
   std::auto_ptr<Column> c = table->add_column("col1", theta::typeHisto);
   Histogram h(1000, 0.0, 1.0);
   for(size_t i=0; i<=1001; ++i){
     h.set(i, i+1);
   }
   //about 8kbytes per row. For 1 mio rows => 8 GB.
   const size_t N = 1000 * 1000;
   BOOST_TEST_CHECKPOINT("about to fill table");
   for(size_t i=0; i<N; ++i){
      if(i % 10000 == 0){
	 BOOST_TEST_MESSAGE("at row " << i << " of " << N);
      }
      table->set_column(*c, h);
      table->add_row();
   }
   BOOST_TEST_CHECKPOINT("deleting database object");
   db.reset();
   BOOST_TEST_CHECKPOINT("database object deleted");
   // for this test to pass, it is enough to come here without any errors ...
   BOOST_REQUIRE(true);
}

BOOST_AUTO_TEST_SUITE_END()
