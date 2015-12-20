#!/bin/sh 
# unset DISPLAY
export CIAO_MAJOR_VER=-4.7

# Use -o option to overrride here. Otherwise, this script fails with an error if a CIAO
# environment is already set up and make stops.
# This is true, evne if the previous CIAO environmet is the same version.
# Due to the >/dev/null the error message is invisible, so it's better to make sure that
# something is set.
. /nfs/cxc/a1/setup/ciao-setup.sh -o > /dev/null
"$@"
