cmake_minimum_required(VERSION 3.15)
project(DemoC C)

set(CMAKE_C_STANDARD 99)

add_executable(BusyMcs Busy/BusyMcs.c)
add_executable(PpMcs Population/PpMcs.c)
