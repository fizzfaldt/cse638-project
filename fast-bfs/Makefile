CXX=cilk++

a.out: bfs.cc bfs.h Makefile fakecilk.h
	$(CXX) -Werror -Wall -O5 -g3 -ggdb -Wextra $< -I /usr/include/x86_64-linux-gnu -B /usr/lib/x86_64-linux-gnu -lrt
