# Note: 
# - all .cpp files in this directory are built as "core-pluygins.so"
# - the .cxx files might or might not be included, depending on the settings

file(GLOB plugin_sources *.cpp)
set(plugin_link_libraries "")

if(psql)
  list(APPEND plugin_sources postgresql_database.cxx)
#  target_link_libraries(plugin pq)
  set(plugin_link_libraries pq)
  include_directories("/usr/include/postgresql/")
endif()

add_library(plugin MODULE ${plugin_sources})
target_link_libraries(plugin libtheta sqlite3 ${plugin_link_libraries} ${Boost_LIBRARIES})

add_custom_command(TARGET plugin POST_BUILD COMMAND mv ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libplugin.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/core-plugins.so)

if(optiontest)
   #build again. For testing, do not build as module to allow testing
   # code to link against it:
   add_library(pluginfortest STATIC ${plugin_sources})
   target_link_libraries(pluginfortest sqlite3)
endif()


