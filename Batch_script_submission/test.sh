#!/bin/bash

echo $0
basename $0
dirname $0
dirname $(readlink -f $0)
SCRIPT_DIR=$(dirname $(readlink -f $0))
echo $SCRIPT_DIR
dirname $SCRIPT_DIR
