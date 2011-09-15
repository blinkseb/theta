add_executable(theta theta.cpp)
target_link_libraries(theta libtheta ${Boost_LIBRARIES})

add_executable(merge merge.cpp)
target_link_libraries(merge sqlite3 libtheta ${Boost_LIBRARIES})
