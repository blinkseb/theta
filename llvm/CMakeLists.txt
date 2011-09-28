execute_process(COMMAND llvm-config --cppflags OUTPUT_VARIABLE llvm_cflags)
execute_process(COMMAND llvm-config --libfiles nativecodegen engine scalaropts ipo OUTPUT_VARIABLE llvm_libs)
string(STRIP "${llvm_cflags}" llvm_cflags)
string(STRIP "${llvm_libs}" llvm_libs)

file(GLOB llvm_sources *.cpp)
add_library(llvm-plugins MODULE ${llvm_sources})
target_link_libraries(llvm-plugins ${Boost_LIBRARIES} libtheta)
add_custom_command(TARGET llvm-plugins POST_BUILD COMMAND mv ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libllvm-plugins.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/llvm-plugins.so)
set_property(TARGET llvm-plugins APPEND PROPERTY LINK_FLAGS "${llvm_libs}")
set_property(TARGET llvm-plugins APPEND PROPERTY COMPILE_FLAGS "${llvm_cflags}")

