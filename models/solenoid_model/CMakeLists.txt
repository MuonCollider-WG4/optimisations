SET(MODULE_NAME fieldtools)
SET(PYTHON_SRCS DerivativesSolenoid.cc)

add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
target_link_libraries( ${MODULE_NAME}
    ${Python3_LIBRARIES}
    m
    z
    "-rdynamic"
)
target_include_directories(${MODULE_NAME} PUBLIC ${Python3_INCLUDE_DIRS})
SET(MY_INSTALL_LOCATION ${PYOPAL_INSTALL_DIR}/pyopal/${INSTALL_DIRECTORY})
MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
#install (TARGETS ${MODULE_NAME} DESTINATION "${MY_INSTALL_LOCATION}")

