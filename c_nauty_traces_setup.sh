#!/bin/bash

cd nauty27r4
./configure --enable-tls
make
git update-index --assume-unchanged nauty27r4/nauty.h
git update-index --assume-unchanged nauty27r4/makefile
git update-index --assume-unchanged nauty27r4/gtools.h
