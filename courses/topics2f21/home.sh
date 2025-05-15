#!/bin/bash
num=$1
echo "Updating homework $num"
cp "$HOME/repos/m401s20/homework/s$num/p$num.pdf" "./hw$num.pdf"
cp "$HOME/repos/m401s20/homework/s$num/s$num.pdf" "./hwsol$num.pdf"

