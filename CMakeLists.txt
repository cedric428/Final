cmake_minimum_required(VERSION 3.9)
project(Final)

set(CMAKE_CXX_STANDARD 11)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(Final main_10.cpp
        BinomialPricer.cpp BinomialPricer.h
        FiniteDifferencePricer.h FiniteDifferencePricer.cpp
        BSpricer.h BSpricer.cpp
        )

target_link_libraries (Final Eigen3::Eigen)
