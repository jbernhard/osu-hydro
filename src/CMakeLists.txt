add_executable(${PROJECT_NAME}
  arsenal.for
  cornelius2.f90
  eos.for
  initialization.for
  input_function.for
  phy_boundary.for
  util.for
  vishnew.for
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
