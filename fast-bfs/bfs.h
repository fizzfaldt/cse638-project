#pragma once
#include <vector>
#include <istream>
#include <iostream>
#include <fstream>
#include <climits>
#include <pthread.h>
#include <inttypes.h>

#define CXX_CHECKS 0
#define SINGLE_THREAD_PARANOID 0
#define CILK_VERIFY 0
#define DEBUG_PRINT 0
#define PARANOID 0
#define QUIT_EARLY 0

#if defined(__cilkplusplus)
// cilk arts
#define _Cilk_for cilk_for
#define _Cilk_sync cilk_sync
#elif defined(__cilk)
// cilkplus
#define cilk_main main
#else
#include "fakecilk.h"
#endif

static const int INFINITY = INT_MAX;
static const int INVALID = INT_MAX;

class Graph;

class Graph {
    private:
        int n;
        int m;
        /* Adjacency list */
        std::vector< std::vector<int> > adj;
        std::vector<int> d;
        std::vector<int> owner;
        int p;
        std::vector< std::vector<int> > q;
        bool opt_c;

    public:
        Graph(bool opt_c);
        void init(int n, int m, std::ifstream &ifs);
        unsigned long long computeChecksum(void);
        int serial_bfs(int s);
        int parallel_bfs(int s);
        int weight(int u, int name);

        static void get_even_split_size_and_offset(int p, int i, int total, int &size, int &offset);

        template<class T>
        static void destructive_serial_prefix_sum(std::vector<T> &v);

        template<class T>
        static void destructive_parallel_prefix_sum(std::vector<T> &v);
    private:
        int find_max_d(void) const;

        template<class T>
        static T destructive_parallel_prefix_sum_up(std::vector<T> &v, size_t start, size_t limit);

        template<class T>
        static void destructive_parallel_prefix_sum_down(std::vector<T> &v, size_t start, size_t limit, bool rightmost_excluded, T partial_sum);

        static int find_index_in_prefix_sum(int value, std::vector<int> &v);
};


class Problem {
    private:
        int n;
        int m;
        int r;
        Graph g;
        std::vector<int> sources;
    public:
        Problem(bool opt_c);
        void init(std::string filename);
        void run(bool parallel);
};

