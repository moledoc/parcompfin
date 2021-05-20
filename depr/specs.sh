#!/bin/sh
# computer specs
neofetch|sed 's/\x1B\[[0-9;]*m//g' > specs.txt
