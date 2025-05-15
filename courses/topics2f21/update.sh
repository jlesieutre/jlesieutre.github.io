#!/bin/bash
num=$1
echo "Updating notes $num"
cp "$HOME/repos/m401s20/notes/s$num/p$num.pdf" "./ws$num.pdf"
cp "$HOME/repos/m401s20/notes/s$num/s$num.pdf" "./wssol$num.pdf"
cp "$HOME/repos/m401s20/notes/s$num/live.pdf" "./slide$num.pdf"

