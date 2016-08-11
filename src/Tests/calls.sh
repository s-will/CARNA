#!/bin/bash

## Test carna

cd Tests

topdir=../$srcdir/.. # top level directory of the locarna source
bin=..

function calltest {
    echo "============================================================"
    echo CALL $*
    echo
    if $* ; then
      echo "======================================== OK"
    else
      echo "======================================== FAILED"
      exit -1
    fi
}

## ========================================
## test carna
##

calltest $bin/carna $topdir/Examples/mouse.fa $topdir/Examples/human.fa -p 0.05 -D 30 -d 60 --time-limit=10000
