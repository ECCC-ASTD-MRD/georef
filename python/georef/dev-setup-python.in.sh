export ECCI_DATA_DIR=/home/sici000/ci_data

source /home/sici000/ci-env/20240722/ppp5/inteloneapi-2024.2.0.sh

# DEV TODO: Change for an official path of librmn with python bindings
. ssmuse-sh -x /home/phc001/Repositories/gitlab.science.gc.ca/RPN-SI/librmn/localinstall
# . ssmuse-sh -d ~phc001/site5/ssm/rmn-with-python.ssmdomain

# Add the build version of the pacakge to PYTHONPATH
export PYTHONPATH=@CMAKE_BINARY_DIR@/python:${PYTHONPATH}


