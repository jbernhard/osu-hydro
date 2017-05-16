add_executable(${PROJECT_NAME}
  arsenal.for
  cornelius2.f90
  eos.for
  initialization.for
  input_function.for
  osu_hydro.for
  phy_boundary.for
  util.for
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
