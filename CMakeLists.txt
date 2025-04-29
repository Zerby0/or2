cmake_minimum_required(VERSION 3.20)

project(tsp)

# genereate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# find sources
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS "src/*.c" "src/*.h" "concorde/*.c" "concorde/*.h")
include_directories(concorde)
set(CMAKE_C_STANDARD 17)

# enable sanitizers in debug mode
option(SANITIZE "Enable sanitizaers in debug builds" ON)
if(SANITIZE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
	add_compile_options(-fomit-frame-pointer -fsanitize=address,undefined)
	add_link_options(-fomit-frame-pointer -fsanitize=address,undefined)
endif()

# dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIRS})

# binary targets
add_executable(tsp ${SRC})
target_link_libraries(tsp PRIVATE ${CPLEX_LIBRARIES})
