add_executable(${PROJECT_NAME}
  arsenal.for
  cornelius2.f90
  eos.for
  initialization.for
  input_function.for
  oscar_output.for
  phy_boundary.for
  vishnew.for
)

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
