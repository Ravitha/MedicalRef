# Script to generate the html information from the coverage tests
# It requires to compile the code with the coverage option (make COV=1) and the
# folder 'outcov' must be created before running the script
lcov --capture --directory . --output-file coverage.info
lcov --remove coverage.info '/usr/include/*' -o coverage_filtered.info
genhtml coverage_filtered.info --output-directory outcov
