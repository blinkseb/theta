execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE root_cflags)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE root_libs)
string(STRIP "${root_cflags}" root_cflags)
string(STRIP "${root_libs}" root_libs)

string(REPLACE "-lNet" "" root_libs "${root_libs}")
string(REPLACE "-lTree" "-lTree -lNet" root_libs "${root_libs}")
#string(REPLACE "-lCore" "" root_libs "${root_libs}")
#set(root_libs "${root_libs} -lCore")



if(root_cflags)
   file(GLOB root_sources *.cpp)
   add_library(root MODULE ${root_sources})
   target_link_libraries(root ${Boost_LIBRARIES} libtheta ${root_libs} -lMinuit2)
   add_custom_command(TARGET root POST_BUILD COMMAND mv ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libroot.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/root.so)
   set_property(TARGET root APPEND PROPERTY COMPILE_FLAGS "${root_cflags}")

   #add_executable(histos histos.cxx)
   #target_link_libraries(histos libtheta sqlite3 ${root_libs} -lMinuit2)
   #set_property(TARGET histos APPEND PROPERTY COMPILE_FLAGS "${root_cflags}")

else(root_cflags)
   message(STATUS "ROOT not found. Not building ROOT dependent stuff")
endif(root_cflags)
