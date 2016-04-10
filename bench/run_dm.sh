#/bin/bash
# usage : find Matrices/alex -name "*.sms" -print0 | xargs -0 -n1  ./run_dm.sh | grep "Dulmage"

printf "%s --- " $1
zcat $1 | ./dm --tabulated
