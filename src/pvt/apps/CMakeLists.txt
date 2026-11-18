add_executable(chflow-pvt pvtApp.cpp)
target_link_libraries(chflow-pvt PRIVATE cf_pvt)

add_executable(chflow-separator separatorApp.cpp)
target_link_libraries(chflow-separator PRIVATE cf_pvt)
