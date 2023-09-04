#!/bin/env bash

# if already sourced, don't do it again
if [ -z "${SOURCED_OPTICS}" ]; then
    echo "Setting up optics"
else
    echo "optics already set up, nothing to do"
    return 0
fi

export SOURCED_OPTICS=1


# add env.sh directory to pythonpath
my_pwd=`pwd`
my_path="${my_pwd}/${BASH_SOURCE[0]}"
my_path=$(dirname "$my_path")
echo "Adding ${my_path} to PYTHONPATH"
export PYTHONPATH=$my_path:$PYTHONPATH
my_python_env="${HOME}/Software/env.sh"

# source software environment
if [ -e  ${my_python_env} ]; then
    source ${my_python_env}
fi