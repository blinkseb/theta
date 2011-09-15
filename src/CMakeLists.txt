file(GLOB libtheta_sources *.cpp)

if(crlibm)
   set(libtheta_sources ${libtheta_sources} ../crlibm/log.c)
endif()

if(generic)
   set(libtheta_sources ${libtheta_sources} log2_dot.cxx)
else()
   set(libtheta_sources ${libtheta_sources} log2_dot.s)
endif()

add_library(libtheta SHARED ${libtheta_sources})
target_link_libraries(libtheta m dl config ${Boost_LIBRARIES})
set_property(TARGET libtheta PROPERTY OUTPUT_NAME theta)

