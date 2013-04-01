#!/bin/bash

#$ -V
#$ -cwd
#$ -q development
#$ -pe 12way 12
#$ -N hw2pa.outb
#$ -o output-$JOB_NAME
#$ -e error-$JOB_NAME
#$ -M rathakrishnanarun@gmail.com
#$ -l h_rt=01:00:00

export PATH=$PATH:$HOME/cilk/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/cilk/lib

export CILK_NPROC=12

#./a.out  /work/01905/rezaul/CSE638/HW1/samples/sample-01-in.txt -a b >sample-01bfs-out_parallel.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/samples/sample-01-in.txt 1 1 0 >sample-01bfs-out_serial.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/samples/sample-02-in.txt 1 1 0 >sample-02bfs-out_parallel.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/samples/sample-02-in.txt -a b >sample-02bfs-out_serial.txt
#./a.out  /work/01945/jesmin/hw2/samples/test-01-in.txt 1 1 0 >test-01bfs-out.txt
#./a.out  /work/01945/jesmin/hw2/samples/test-02-in.txt 1 1 0 >test-02bfs-out.txt
#./a.out  /work/01945/jesmin/hw2/samples/test-03-in.txt 1 1 0 >test-03bfs-out.txt
#./a.out  ../graph-10M-100M-1k.in 1 1 0 >graph-10M-100M-1k-bfsMIT-out.txt
#./a.out  /work/01945/jesmin/hw2/samples/BFSFinalVersions/graph-10M-1B-1k.in 1 1 0 >graph-10M-1B-1k.txt

#cilkview -trials all 12 ./a.out  /work/01945/jesmin/hw2/samples/kkt_power.txt 1 1 0 > kkt_powerpbfsbT12-out.txt
#cilkview -trials all 12 ./a.out  /work/01945/jesmin/hw2/samples/Freescale1.txt 1 1 0 > Freescale1bfsbT12-out.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/kkt_power-in.txt 1 1 0 > kkt_powerpbfsbT12-out.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/cage14-in.txt 1 1 0 > cage14-out.txt 2> cage14-err.txt
#./a.out /work/01905/rezaul/CSE638/HW1/turn-in/cage15-in.txt 1 1 0 > cage15-out.txt 2> cage15-err.txt
./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/wikipedia-in.txt 1 1 0 > wikipedia-out.txt 2>wikipedia-err.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/freescale-in.txt 1 1 0 > freescale-out.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/rmat100M-in.txt 1 1 0 > rmat100M-out.txt
#./a.out  /work/01905/rezaul/CSE638/HW1/turn-in/rmat1B-in.txt 1 1 0 > rmat1B-out.txt
