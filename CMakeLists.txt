cmake_minimum_required(VERSION 2.6)

project("DHLLDV-hpc")

add_definitions(-std=c++11)
include_directories(/usr/local/include/)
link_directories(/usr/local/lib/)

set(BOOST_INCLUDEDIR ../boost_1_64_0/)

find_package(Boost REQUIRED)

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(test_LDV test_LDV.cpp
                        waterData.cpp 
                        limitDepositVelocity.cpp 
                        particleUtils.cpp 
                        fluidutils.cpp )
                        
add_executable(test_DHLLDV_class test_DHLLDV_class.cpp
                                DHLLDV.cpp
                                waterData.cpp 
                                particleUtils.cpp 
                                fluidutils.cpp
                                core/utils.cpp )
                                
add_executable(test_slipRatio test_slipRatio.cpp
                                waterData.cpp
                                slipRatio.cpp 
                                particleUtils.cpp 
                                fluidutils.cpp 
                                limitDepositVelocity.cpp )
                                
add_executable(imitateXlsx imitateXlsx.cpp 
                            DHLLDV.cpp
                            particleUtils.cpp
                            fluidutils.cpp
                            core/utils.cpp )
                            
target_link_libraries(imitateXlsx xlsxwriter)