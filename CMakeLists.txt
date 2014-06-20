project(breast-cancer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")
add_executable(breast-cancer
    main.cpp
    samples.cpp
    p_statistics.cpp
    classificator.cpp
    )
