#!/bin/bash
sed -e 's/\([^\t\n\s:]*\):[^\t\n\s]*[\t\n\s]/\1\t/g' $1 | sed -e 's/\([^\t\n\s:]*\):[^\t\n\s]*$/\1/g' 
