#/bin/bash
# usage : find Matrices/alex -name "*.sms" -print0 | xargs -0 -n1  ./run_dm.sh | grep "Dulmage"

OUTPUT_PATH=../result/
printf "%s --- " $1
zcat $1 | ./../../gbla/tools/dump_matrix -n -s - | ./demo_cheap > $OUTPUT_PATH/gbla.sms
