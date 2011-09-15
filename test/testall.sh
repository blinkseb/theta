#!/usr/bin/env bash

#run this from the theta main directory!

[ -d test ] || { echo "not in theta root directory!"; exit 1; }

. test/lib.sh

[ -x root/create_testhistos ] && root/create_testhistos
execute_checked bin/test
rm -f testhistos.root

fail=0

for i in test/test-stat/*.py; do
   i=`basename $i`
   echo "START `date +%s`;     `date -R`     $i"
   export output=$( cd test/test-stat; ./$i 2>&1 )
   echo "output: ${output}"
   output_nonpass=$( echo "$output" | grep -v "PASS:" )
   if [ -n "$output_nonpass" ]; then
       echo -e "\n\e[0;31m $i FAILED: ${output_nonpass}\e[m \n";
       fail=$(($fail+1))
   fi
   echo "$output" | awk "{print \"[$i]\", \$_}"
   echo "END `date +%s`;     `date -R`     $i"
done

echo "Failures: ${fail}"
[ $fail -eq 0 ] || exit 1

