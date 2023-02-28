#!/bin/bash
source ~/.bashrc
nohup ./job.sh > nohup.out 2>&1 &
echo $!
