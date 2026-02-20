cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(UserActions)

# Load ActRoot custom script
set(ACTROOT $ENV{ACTROOT})
if(NOT ACTROOT)
    message(FATAL_ERROR "ACTROOT environment variable is not set")
endif()

include(${ACTROOT}/Scripts.cmake)

# First user action
add_userlibrary(
    NAME RecRANSAC
    SOURCES RecRANSAC.cxx   # solo el .cxx
    LINK ActAlgorithm        # ActAlgorithm debe estar en $ACTROOT/install/lib
)
add_userlibrary(NAME Task0 SOURCES Task0.cxx LINK ActAlgorithm)

