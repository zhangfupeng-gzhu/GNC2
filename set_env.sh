#!/bin/bash
export GWPATH=$(pwd)
export PATH=$PATH:$GWPATH"/source/main/"
export PYTHONPATH=$GWPATH"/pycom":$PYTHONPATH
ulimit -s 2048000