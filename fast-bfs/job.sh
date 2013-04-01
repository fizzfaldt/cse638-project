#!/bin/bash

#$ -V
#$ -cwd
#$ -q development
#$ -pe 12way 12
#$ -N hw2pbfsb
#$ -o output-$JOB_NAME
#$ -e error-$JOB_NAME
#$ -M rathakrishnanarun@gmail.com
#$ -l h_rt=01:00:00

export PATH=$PATH:$HOME/cilk/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/cilk/lib

export CILK_NPROC=12

#./bfs -f /work/01905/rezaul/CSE638/HW1/samples/sample-01-in.txt -a b >sample-01bfs-out_parallel.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/samples/sample-01-in.txt -a p >sample-01bfs-out_serial.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/samples/sample-02-in.txt -a p >sample-02bfs-out_parallel.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/samples/sample-02-in.txt -a b >sample-02bfs-out_serial.txt
#./bfs -f /work/01945/jesmin/hw2/samples/test-01-in.txt -a p >test-01bfs-out.txt
#./bfs -f /work/01945/jesmin/hw2/samples/test-02-in.txt -a p >test-02bfs-out.txt
#./bfs -f /work/01945/jesmin/hw2/samples/test-03-in.txt -a p >test-03bfs-out.txt
#./bfs -f ../graph-10M-100M-1k.in -a p >graph-10M-100M-1k-bfsMIT-out.txt
#./bfs -f /work/01945/jesmin/hw2/samples/BFSFinalVersions/graph-10M-1B-1k.in -a p >graph-10M-1B-1k.txt

#cilkview -trials all 12 ./bfs -f /work/01945/jesmin/hw2/samples/kkt_power.txt -a p> kkt_powerpbfsbT12-out.txt
#cilkview -trials all 12 ./bfs -f /work/01945/jesmin/hw2/samples/Freescale1.txt -a p> Freescale1bfsbT12-out.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/kkt_power-in.txt -a p> kkt_powerpbfsbT12-out.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/cage14-in.txt -a p> cage14-out.txt
./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/cage15-in.txt -a p> cage15-out.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/wikipedia-in.txt -a p> wikipedia-out.txt
./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/freescale-in.txt -a p> freescale-out.txt
./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/rmat100M-in.txt -a p> rmat100M-out.txt
#./bfs -f /work/01905/rezaul/CSE638/HW1/turn-in/rmat1B-in.txt -a p> rmat1B-out.txt
