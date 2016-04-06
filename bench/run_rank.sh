#/bin/bash
# usage : find Matrices/alex -name "*.sms" -print0 | xargs -0 -n1  ./run_dm.sh | grep "Dulmage"

printf "%s --- " $1
if gzcat $1 | ./rank --max-time 10 &> /dev/null; then
    printf "OK\n"
else
    printf "FAIL\n"
fi
#convert plop.ppm $1_dm.png

