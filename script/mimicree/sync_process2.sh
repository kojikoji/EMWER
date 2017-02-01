#!/bin/bash
sed -e 's/\([^\t\n:]*\):[^\t\n]*[\t\n]/\1\t/g' $1 | sed -e 's/\([^\t\n:]*\):[^\t\n]*$/\1/g' 
