rm gnu-baseline.txt
touch gnu-baseline.txt
echo "******** Start Base Line ********"
for i in `seq 1 10`
do
    echo $i
   ./GNU-baseline/tsvc_novec_default >>gnu-baseline.txt
done
rm gnu-sse.txt
touch gnu-sse.txt
echo "******** Start SSE ********"
for i in `seq 1 5`
do
    echo $i
   ./GNU-sse/tsvc_novec_default >>gnu-sse.txt
done
rm gnu-avx.txt
touch gnu-avx.txt
echo "******** Start AVX ********"
for i in `seq 1 10`
do
    echo $i
   ./GNU-avx/tsvc_novec_default >>gnu-avx.txt
done
