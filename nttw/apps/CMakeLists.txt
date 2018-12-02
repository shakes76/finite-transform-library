# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (machine_properties)
add_subdirectory (ntt)
add_subdirectory (intt)
add_subdirectory (operate)
add_subdirectory (farey_points)

if(BUILD_ITK_SUPPORT)
  add_subdirectory (pgm_convert)
endif(BUILD_ITK_SUPPORT)