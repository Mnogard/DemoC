cmake_minimum_required(VERSION 3.15)
project(DemoC C)

set(CMAKE_C_STANDARD 99)

add_executable(BusyMcs Busy/BusyMcs.c)
add_executable(BusyHot Busy/BusyHot.c)
add_executable(BusyQu Busy/BusyQu.c)
add_executable(BusyBt Busy/BusyBt.c)
add_executable(DemoBusyMcs Busy/DemoBusyMcs.c)


add_executable(PpMcs Population/PpMcs.c)
add_executable(PpHot Population/PpHot.c)
add_executable(Hope Hope/Hope.c)

add_executable(MppMcs MemoryPp/MppMcs.c)

add_executable(Demo Demo.c)