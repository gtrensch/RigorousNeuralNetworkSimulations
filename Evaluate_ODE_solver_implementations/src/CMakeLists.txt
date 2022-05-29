cmake_minimum_required (VERSION 3.9)

# + + + set a project name, it is referenced in other cmake function calls
project( izhikevich )

# + + + set cache variables for a version number
set( VERSION_MAJOR 1 )
set( VERSION_MINOR 0 )

# + + + add executable target to the project:   add_executable( <target_name> source1 [source 2] ...)
add_executable( izhikevich izhikevich izhikevich.cpp )

# + + + specify  include directories
include_directories( izhikevich /usr/include/plplot )

# + + + ... and libraries
target_link_libraries( izhikevich libplplot.so libplplotcxx.so gsl gslcblas )
