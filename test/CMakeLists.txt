
file(GLOB test_sources test*.cpp utils.cpp)
add_executable(test-exe ${test_sources})
target_link_libraries(test-exe libtheta)

add_library(libtest MODULE test-ex2.cxx  test-ex3.cxx  test-ex.cxx)

file(GLOB plugintest_sources plugintest*.cpp utils.cpp test.cpp)
add_executable(plugintest-exe ${plugintest_sources})
target_link_libraries(plugintest-exe libtheta pluginfortest)

