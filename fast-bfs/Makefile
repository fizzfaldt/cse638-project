CXX=cilk++
DEBUG=0
ifeq ($(DEBUG),1)
    OPTIMIZE=-O0
else
    OPTIMIZE=-O5
endif

a.out: bfs.cc bfs.h Makefile fakecilk.h
	$(CXX) -Werror -Wall $(OPTIMIZE) -g3 -ggdb -Wextra $< -I /usr/include/x86_64-linux-gnu -B /usr/lib/x86_64-linux-gnu -lrt
