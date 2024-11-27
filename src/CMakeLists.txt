add_library(Precision precision.f95)

# Define configuration
add_library(Config config.f95)
add_dependencies(Config Precision)


add_library(TRI tridiagonal.f95)
add_dependencies(TRI Precision)
target_link_libraries(TRI PUBLIC Precision)

add_library(MD multidim.f95)
add_dependencies(MD Precision TRI)
target_link_libraries(MD PUBLIC Precision TRI)




# Make final executable
add_executable(multidim main.f95)

target_link_libraries(multidim PUBLIC 
	Config Precision TRI MD)


# Move executable to thr working directiory
install(TARGETS multidim DESTINATION ${CMAKE_BIN_DIR})
