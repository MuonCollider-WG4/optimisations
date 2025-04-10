cmake_minimum_required(VERSION 3.10) # probably
project(ProjectName)

IF (DEFINED ENV{SOLENOID_OPTICS_SRC_DIR})
	MESSAGE(STATUS "USING $ENV{SOLENOID_OPTICS_DIR}")
ELSE()
	MESSAGE(FATAL_ERROR "SOLENOID_OPTICS_SRC_DIR was not set - check 'source env.sh'")
	RETURN()
ENDIF()

SET(MODULE_NAME pulsed_rf)
SET(PYTHON_SRCS PulsedRF.cc $ENV{SOLENOID_OPTICS_SRC_DIR}/models/solenoid_model/DerivativesSolenoid.cc PyPulsedRF.cc)

add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
include_directories(${Python_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
target_include_directories(${MODULE_NAME} PUBLIC ${Python_INCLUDE_DIR} $ENV{SOLENOID_OPTICS_INCLUDE_DIR})
target_link_libraries( ${MODULE_NAME}
    ${Python_LIBRARY}
    ${GSL_LIBRARY}
    m
    z
    "-rdynamic"
    "-std=c++2a"
)
SET(MY_INSTALL_LOCATION ${PYOPAL_INSTALL_DIR}/pyopal/${INSTALL_DIRECTORY})
MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
#install (TARGETS ${MODULE_NAME} DESTINATION "${MY_INSTALL_LOCATION}")

