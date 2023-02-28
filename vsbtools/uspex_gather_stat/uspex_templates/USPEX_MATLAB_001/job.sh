#!/bin/sh

while [ ! -f ./USPEX_IS_DONE ]; do
date >> log
matlab < USPEX.m >> log -singleCompThread
sleep 2
done
