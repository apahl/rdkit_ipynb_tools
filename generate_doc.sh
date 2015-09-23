#!/usr/bin/env bash
pdoc --html --overwrite ./

CWD=${PWD##*/}
rm -rf doc/
mv $CWD doc