cmake_minimum_required(VERSION 3.20)
project(FEMVisualizer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Qt6 COMPONENTS Core Widgets OpenGL OpenGLWidgets REQUIRED)
find_package(GTest CONFIG REQUIRED)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)