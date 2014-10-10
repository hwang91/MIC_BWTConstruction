
./partition ~/Documents/Data/read1.2M.fa 100 100

 time for i in {0..255}; do (./a.out BWT_$i partitionCount.ini $i); done
