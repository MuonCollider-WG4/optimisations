SET(MODULE_NAME derivatives_solenoid)
SET(PYTHON_SRCS DerivativesSolenoid.cc PyDerivativesSolenoid.cc)

add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
include_directories(${Python_INCLUDE_DIR})
target_include_directories(${MODULE_NAME} PUBLIC ${Python_INCLUDE_DIR})
target_link_libraries( ${MODULE_NAME}
    ${Python_LIBRARY}
    m
    z
    "-rdynamic"
    "-std=c++2a"
)
SET(MY_INSTALL_LOCATION ${PYOPAL_INSTALL_DIR}/pyopal/${INSTALL_DIRECTORY})
MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
#install (TARGETS ${MODULE_NAME} DESTINATION "${MY_INSTALL_LOCATION}")

