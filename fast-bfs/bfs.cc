/*
 * n m r
 * //m
 * //r (sources, repeat algo for this)
 * m lines: u v (sorted by u) <- directed edge
 * r lines sources: s (
 */
#include <assert.h>
#include <queue>
#include <utility>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#ifdef __cilkplusplus
#include <cilk.h>
#include <reducer_opadd.h>
#include <reducer_opand.h>
#include <reducer_max.h>
//#include <cilk_mutex.h>
#endif

#include "bfs.h"

using namespace std;
cilk::context ctx;

Graph::Graph(bool opt_c) :
    n(0),
    m(0),
    adj(),
    d(),
    owner(),
    p(ctx.get_worker_count()),
    opt_c(opt_c) {
    srandom(time(NULL));
}

void Graph::init(int n, int m, ifstream &ifs) {
    adj.clear();
    adj.resize(n);
    d.resize(n);
    owner.resize(n);
    q.resize(p);

    this->n = n;
    this->m = m;

    for (int i = 0; i < m; ++i) {
        int u, v;
        ifs >> u >> v;
        assert(u>0);
        assert(v>0);
        adj[u-1].push_back(v-1);
    }
}

#define FOR_EACH_GAMMA(v, vec) for (vector<int>::iterator v = (vec).begin(); v != (vec).end(); ++(v))

unsigned long long Graph::computeChecksum(void) {
    cilk::reducer_opadd<unsigned long long> chksum;

    cilk_for (int i = 0; i < n; ++i) {
        chksum += d[i] == INFINITY ? n : d[i];
    }
#if CILK_VERIFY 
    int mind = INT_MAX;
    int maxd = INT_MIN;
    unsigned long long chksum_ser = 0;
    for (int i = 0; i < n; ++i) {
        assert(d[i] >= 0);
        int change = d[i] == INFINITY ? n : d[i];
        chksum_ser += change;
        mind = std::min(mind, change);
        maxd = std::max(maxd, change);
    }
    fprintf(stderr, "min=%d, max=%d\n", mind, maxd);
    fflush(stderr);

    assert(chksum.get_value() == chksum_ser);
#endif

    return chksum.get_value();
}

int Graph::serial_bfs(int s) {
    for (int u = 0; u < n; ++u) {
        d[u] = INFINITY;
    }
    d[s] = 0;
    queue<int> Q;
    Q.push(s);
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        FOR_EACH_GAMMA(v, adj[u]) {
            if (d[*v] == INFINITY) {
                d[*v] = d[u] + 1;
                Q.push(*v);
            }
        }
    }
    int maxd = 0;
    for (int u = 0; u < n; ++u) {
    	if(d[u] == INFINITY) continue;
        maxd = max(maxd, d[u]);
    }
    return maxd;
}

void Graph::get_even_split_size_and_offset(int i, int total, int *size, int *offset) const {
    if (i < total % p) {
        *size = (total + p - 1) / p;  //ceiling
    } else {
        *size = total / p;  // floor
    }
    *offset = i * (total / p) + min(i, total % p);
}

int Graph::find_index_in_prefix_sum(int value, vector<int> &v) const {
    int low = 0;
    int high = v.size() - 1;
    assert(high > low);
    assert(v.back() >= value);

    while (high > low) {
        int mid = (high + low) / 2;
        if (v[mid] < value) {
            // Search in upper half
            low = mid + 1;
            assert(low <= high);
        } else if (v[mid] > value && mid > low && v[mid-1] >= value) {
            // Search in lower half
            high = mid - 1;
        } else {
            // Found it.  Break out for single return statement.
            low = high = mid;
        }
    }
    assert(high == low);
    assert(v[low] >= value);
    return low;
}

void Graph::destructive_serial_prefix_sum(std::vector<int> &v) const {
    int sum = 0;
    for (int i = 0; i < (int)v.size(); i++) {
        sum += v[i];
        v[i] = sum;
    }
}

int Graph::destructive_parallel_prefix_sum_up(std::vector<int> &v, int start, int limit) const {
    const int size = limit - start;
    assert(size > 0);
    if (size == 1) {
        return v[start];
    }
    int x = cilk_spawn destructive_parallel_prefix_sum_up(v, start, (start+limit)/2);
    int y = destructive_parallel_prefix_sum_up(v, (start+limit)/2, limit);
    return (v.back() = x+y);
}

void Graph::destructive_parallel_prefix_sum_down(std::vector<int> &v, int start, int limit, int partial_sum) const {
    const int size = limit - start;
    assert(size > 0);
    if (size == 1) {
        v[start] += partial_sum;
        return;
    }
    int sum_left = v[(start+limit)/2 - 1];
    cilk_spawn destructive_parallel_prefix_sum_down(v, start, (start+limit)/2, partial_sum);
    destructive_parallel_prefix_sum_down(v, (start+limit)/2, limit, partial_sum + sum_left);
}

void Graph::destructive_parallel_prefix_sum(std::vector<int> &v) const {
    if (v.size() > 0) {
        destructive_parallel_prefix_sum_up(v, 0, v.size());
        destructive_parallel_prefix_sum_down(v, 0, v.size(), 0);
        // Optional:  cilk_sync can be moved to last line of destructive_parallel_prefix_sum_down
        cilk_sync;
    }
}

int Graph::parallel_bfs(int s) {
    cilk_for(int u = 0; u < n; ++u) {
        d[u] = INFINITY;
        if (opt_c) {
            owner[u] = INVALID;
        }
    }
#if DEBUG_PRINT
    fprintf(stderr, "setting d[%d] source = 0. n=%d\n", s, n);
    fflush(stderr);
#endif
    d[s] = 0;
    if (opt_c) {
        owner[s] = 0;
    }
    int input_vertexes[n];
    input_vertexes[0] = s;
    int num_input_vertexes = 1;;

    vector<int> prefix_sum_degrees[p];
    vector<int> prefix_sum_degrees_per_core;
    prefix_sum_degrees_per_core.resize(p, 0);

    int num_degrees[p];
    int offset_degrees[p];
    int num_vertexes[p];
    int offset_vertexes[p];
    vector<int> queue_size;
    queue_size.resize(p, 0);
    q.resize(p);

    while (num_input_vertexes > 0) {
        cilk_for (int i = 0; i < p; i++) {
            get_even_split_size_and_offset(i, num_input_vertexes, &num_vertexes[i], &offset_vertexes[i]);

            // Generate prefix sum for out-degrees in core i's portion of input_vertexes
            prefix_sum_degrees[i].clear();
            prefix_sum_degrees[i].reserve(num_vertexes[i]);
            for (int j = 0; j < num_vertexes[i]; j++) {
                const int u = input_vertexes[j+offset_vertexes[i]];
                prefix_sum_degrees[i].push_back(adj[u].size());
            }
            destructive_serial_prefix_sum(prefix_sum_degrees[i]);
            assert((int)prefix_sum_degrees[i].size() == num_vertexes[i]);
            prefix_sum_degrees_per_core[i] = prefix_sum_degrees[i].back();
        }
        destructive_parallel_prefix_sum(prefix_sum_degrees_per_core);
        int total_degrees = prefix_sum_degrees_per_core.back();
        cilk_for (int i = 0; i < p; i++) {
            get_even_split_size_and_offset(i, total_degrees, &num_degrees[i], &offset_degrees[i]);

            const int sublist = find_index_in_prefix_sum(offset_degrees[i], prefix_sum_degrees_per_core);
            const int degrees_in_earlier_lists = (i==0) ? 0 : prefix_sum_degrees_per_core[i-1];
            const int vertex_in_sublist = find_index_in_prefix_sum(offset_degrees[i] - degrees_in_earlier_lists, prefix_sum_degrees[sublist]);
            const int degrees_in_earlier_vertexes = (vertex_in_sublist==0) ? 0 : prefix_sum_degrees[sublist][vertex_in_sublist-1];

            const int starting_degree = offset_degrees[i] - degrees_in_earlier_lists - degrees_in_earlier_vertexes;
            assert(starting_degree >= 0);
            const int starting_vertex = vertex_in_sublist + offset_vertexes[i];

            int remaining = num_degrees[i];
            int vertex = starting_vertex;
            int degree = starting_degree;
            q[i].clear();
            while (remaining > 0) {
                int u = input_vertexes[vertex];
                int limit = min(remaining + degree, (int)adj[u].size());
                for (int j = degree; j < limit; j++) {
                    int v = adj[u][j];
                    if (d[v] == INFINITY) {
                        owner[v] = i;  // Dedup/duplicate expansion optimization
                        d[v] = d[u] + 1;
                        q[i].push_back(v);
                    }
                }
                remaining -= limit - degree;
                degree = 0;
                vertex++;
            }
            assert(remaining == 0);
        }
        // Delete duplicate items from queues.
        cilk_for (int i = 0; i < p; i++) {
            //TODO: This can be balanced better by doing searches (including 'next queue' searches)
            int valid = 0;
            for (int index = 0; index < (int)q[i].size(); index++) {
                const int u = q[i][index];
                if (owner[u] == i) {
                    q[i][valid++] = u;
                }
            }
            q[i].resize(valid);

            //TODO: If we disable (possibly by option) dedupping, we still need the next line.
            //It can be moved up to prev cilk_for if no dedupping is done.
            queue_size[i] = q[i].size();
        }
        destructive_parallel_prefix_sum(queue_size);
        num_input_vertexes = queue_size.back();
        assert(num_input_vertexes <= n);
        cilk_for (int i = 0; i < p; i++) {
            const int offset = (i==0) ? 0 : queue_size[i-1];
            for (int j = 0; j < (int)q[i].size(); j++) {
                input_vertexes[offset+j] = q[i][j];
            }
        }
    }

    cilk::reducer_max<int> maxd(0);
    cilk_for (int u = 0; u < n; ++u) {
    	if(d[u] == INFINITY) continue;
        maxd = cilk::max_of(maxd, d[u]);
    }
#if CILK_VERIFY
    int maxd_ser = 0;
    for (int u = 0; u < n; ++u) {
    	if(d[u] == INFINITY) continue;
        maxd_ser = max(maxd_ser, d[u]);
    }
    assert(maxd_ser == maxd.get_value());
#endif
    return maxd.get_value();
}

Problem::Problem(bool opt_c) :
    n(0),
    m(0),
    r(0),
    g(opt_c),
    sources() {
}

void Problem::init(string filename) {
    ifstream ifs(filename.c_str());
    ifs >> n >> m >> r;

    g.init(n, m, ifs);
    for (int i = 0; i < r; ++i) {
        int s;
        ifs >> s;
        assert(s>0);
        sources.push_back(s-1);
    }
    ifs.close();
}

void Problem::run(bool parallel) {
    for (int s = 0; s < r; ++s) {
        int maxd;
        if (parallel) {
            maxd = g.parallel_bfs(sources[s]);
        } else {
            maxd = g.serial_bfs(sources[s]);
        }
        printf("%d %llu\n", maxd, g.computeChecksum());
#if QUIT_EARLY 
        fprintf(stderr, "QUITTING\n");
        break;
#endif
    }
}

int cilk_main(int argc, char *argv[]) {
    assert(argc == 6);
    struct timespec start_t;
    struct timespec end_t;
    int error;
    error = clock_gettime(CLOCK_REALTIME, &start_t);
    assert(error==0);

    string file_name(argv[1]);
    bool opt_c = atoi(argv[2]);
    bool do_parallel = atoi(argv[5]);

    if (!opt_c) {
        fprintf(stderr, "Disabling option c (dedup) is not yet supported\n");
        exit(1);
    }
    Problem p(opt_c);
    p.init(file_name);
    p.run(do_parallel);

    error = clock_gettime(CLOCK_REALTIME, &end_t);
    assert(error==0);

    double start_d = start_t.tv_sec + start_t.tv_nsec / 1e9;
    double end_d = end_t.tv_sec + end_t.tv_nsec / 1e9;

    fprintf(stderr, "Time taken: %f\n", end_d - start_d);
//    p.init(filename);
    return 0;
}

