#############################################################################
#
# The script sets up the environmental varibales for Foerster Boxes:
# should work for bash, dash, zsh
#
# a) LIMCON_POST_ROOT
# b) PYTHONPATH
#
##############################################################################

# Get the base directory of this script
if test -n "$BASH" ; then script=$BASH_SOURCE
elif test -n "$TMOUT"; then script=${.sh.file}
elif test -n "$ZSH_NAME" ; then script=${(%):-%x}
elif test ${0##*/} = dash; then x=$(lsof -p $$ -Fn0 | tail -1); script=${x#n}
else script=$0
fi

# use base directory of the script as LIMCON_POST_ROOT
LIMCON_POST_ROOT=$(cd `dirname $script` && pwd)/
export LIMCON_POST_ROOT

echo "Set LIMCON_POST_ROOT to: " ${LIMCON_POST_ROOT}

# Set the python path for all subprojects of Foerster magnetometer
export PYTHONPATH=${LIMCON_POST_ROOT}:${LIMCON_POST_ROOT}/lib/:$PYTHONPATH

