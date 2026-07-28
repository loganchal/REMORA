# - Find PNetCDF
# Find the native PNetCDF includes and library
#
#  PNETCDF_INCLUDE_DIRS - where to find pnetcdf.h, etc
#  PNETCDF_FOUND        - True if PNetCDF found
#
# The following are not for general use and are included in
# PNETCDF_LIBRARIES if the corresponding option above is set.
#
#  PNETCDF_LIBRARIES      - only the libraries (without the '-l')
#  PNETCDF_LINK_LIBRARIES - the libraries and their absolute paths
#  PNETCDF_LDFLAGS        - all required linker flags
#
# Normal usage would be:
#  find_package (PNetCDF REQUIRED)
#  target_link_libraries (target_name PUBLIC ${PNETCDF_LINK_LIBRARIES})

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  # Already in cache, be silent
  set (PNETCDF_FIND_QUIETLY TRUE)
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
  pkg_check_modules(PNETCDF IMPORTED_TARGET pnetcdf)
endif ()

# Fall back to a direct search when pkg-config is unavailable or the .pc file
# is not on PKG_CONFIG_PATH. PNETCDF_DIR (or the environment variable of the
# same name) may be set to the installation prefix.
if (NOT PNETCDF_LINK_LIBRARIES)
  find_path (PNETCDF_INCLUDE_DIRS pnetcdf.h
             HINTS ${PNETCDF_DIR} ENV PNETCDF_DIR
             PATH_SUFFIXES include)
  find_library (PNETCDF_LINK_LIBRARIES NAMES pnetcdf
                HINTS ${PNETCDF_DIR} ENV PNETCDF_DIR
                PATH_SUFFIXES lib lib64)
  if (PNETCDF_LINK_LIBRARIES)
    set (PNETCDF_LIBRARIES pnetcdf)
  endif ()
endif ()

# handle the QUIETLY and REQUIRED arguments and set PNETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PNetCDF DEFAULT_MSG PNETCDF_LIBRARIES PNETCDF_LINK_LIBRARIES PNETCDF_INCLUDE_DIRS)

mark_as_advanced (PNETCDF_LIBRARIES PNETCDF_INCLUDES)

