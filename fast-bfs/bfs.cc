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

#if defined(__cilkplusplus)
#include <cilk.h>
#include <reducer_opadd.h>
#include <reducer_opand.h>
#include <reducer_max.h>
#elif defined(__cilk)
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_opand.h>
#include <cilk/reducer_max.h>
#endif

//#include <cilk_mutex.h>

#include "bfs.h"
#define DEBUG_ARRAY 0

using namespace std;
static int get_num_workers() {
#if defined(__cilkplusplus) || !defined(__cilk)
// cilk arts
    static cilk::context ctx;
    return ctx.get_worker_count();
#else
    return __cilkrts_get_nworkers();
#endif
}

Graph::Graph(bool opt_c) :
    n(0),
    m(0),
    adj(),
    d(),
    owner(),
    p(get_num_workers()),
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

    _Cilk_for (int i = 0; i < n; ++i) {
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
int Graph::get_slot_number(int p, int total, int degree){
	//each of the first early_slots will have total/p + 1 elements, total elements in these slots is early_elements
	//if degree > early_elements first early_elements will fit in early_slots each of size total/p +1
	//the remaining_elements will fit in slots of size total/p
	if(p >= total) return degree-1;
	int early_slots = total % p;
	int early_elements = early_slots * ((total / p) + 1);
	int remaining_elements = degree - early_elements;
	if(degree <= early_elements)
		return (degree -1)/((total / p)+1);
	return early_slots + (remaining_elements - 1)/(total /p);
}
void Graph::get_even_split_size_and_offset(int p, int i, int total, int &size, int &offset) {
    if (i < total % p) {
        size = (total + p - 1) / p;  //ceiling
    } else {
        size = total / p;  // floor
    }
    offset = i * (total / p) + min(i, total % p);
}

int Graph::find_index_in_prefix_sum(int value, vector<int> &v) {
    int low = 0;
    int high = v.size() - 1;
    assert(!v.empty());
    assert(v.back() >= value);
    assert(value >= 0);

    while (high > low) {
        int mid = (high + low) / 2;
        if (v[mid] < value) {
            // Search in upper half
            low = mid + 1;
            assert(low <= high);
        } else if (v[mid] >= value && mid > low && v[mid-1] >= value) {
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

template<class T>
void Graph::destructive_serial_prefix_sum(std::vector<T> &v) {
    T sum = 0;
    for (size_t i = 0; i < v.size(); i++) {
        sum += v[i];
        v[i] = sum;
    }
}

int Graph::destructive_parallel_prefix_sum_up(int * v, size_t start, size_t limit) {
    assert(limit > start);
    const size_t size = limit - start;
    if (size == 1) {
        return v[start];
    }
    int x, y;
    x = _Cilk_spawn destructive_parallel_prefix_sum_up(v, start, (start+limit)/2);
    y = destructive_parallel_prefix_sum_up(v, (start+limit)/2, limit);
    _Cilk_sync;
    return (v[limit-1] = x+y);
}

template<class T>
T Graph::destructive_parallel_prefix_sum_up(std::vector<T> &v, size_t start, size_t limit) {
    assert(limit > start);
    const size_t size = limit - start;
    if (size == 1) {
        return v[start];
    }
    T x, y;
    x = _Cilk_spawn destructive_parallel_prefix_sum_up(v, start, (start+limit)/2);
    y = destructive_parallel_prefix_sum_up(v, (start+limit)/2, limit);
    _Cilk_sync;
    return (v[limit-1] = x+y);
}

void Graph::destructive_parallel_prefix_sum_down(int * v, size_t start, size_t limit, bool rightmost_excluded, int partial_sum) {
    assert(limit > start);
    const size_t size = limit - start;
    if (size == 1) {
        v[start] += partial_sum;
        return;
    }
    int sum_left = v[(start+limit)/2 - 1];
    size_t mid = (start+limit)/2;
    // Do entire left side.
    _Cilk_spawn destructive_parallel_prefix_sum_down(v, start, mid, false, partial_sum);
    // Rightmost element just takes partial_sum and is then excluded.
    if (!rightmost_excluded) {
        v[limit-1] += partial_sum;
    }
    assert(limit > mid);
    // Do right half (excluding last element)
    if (limit - mid > 1) {
        destructive_parallel_prefix_sum_down(v, mid, limit, true, partial_sum + sum_left);
    }
    _Cilk_sync;
}


template<class T>
void Graph::destructive_parallel_prefix_sum_down(std::vector<T> &v, size_t start, size_t limit, bool rightmost_excluded, T partial_sum) {
    assert(limit > start);
    const size_t size = limit - start;
    if (size == 1) {
        v[start] += partial_sum;
        return;
    }
    T sum_left = v[(start+limit)/2 - 1];
    size_t mid = (start+limit)/2;
    // Do entire left side.
    _Cilk_spawn destructive_parallel_prefix_sum_down(v, start, mid, false, partial_sum);
    // Rightmost element just takes partial_sum and is then excluded.
    if (!rightmost_excluded) {
        v[limit-1] += partial_sum;
    }
    assert(limit > mid);
    // Do right half (excluding last element)
    if (limit - mid > 1) {
        destructive_parallel_prefix_sum_down(v, mid, limit, true, partial_sum + sum_left);
    }
    _Cilk_sync;
}

__attribute__((__used__))
static void printv(vector<uint64_t> v, char *name) {
    printf("V_%s=\n", name);
    for (size_t i = 0; i < v.size(); i++) {
        printf("\t[%d] = [%lu]\n", (int)i, v[i]);
    }
}

void Graph::destructive_parallel_prefix_sum(int * v, size_t limit) {
    if (limit > 0) {
        destructive_parallel_prefix_sum_up(v, 0, limit);
        destructive_parallel_prefix_sum_down(v, 0, limit, false, 0);
    }
}

template<class T>
void Graph::destructive_parallel_prefix_sum(std::vector<T> &v) {
    if (!v.empty()) {
        destructive_parallel_prefix_sum_up(v, 0, v.size());
        destructive_parallel_prefix_sum_down(v, 0, v.size(), false, static_cast<T>(0));
    }
}

int Graph::parallel_bfs(int s) {
    _Cilk_for(int u = 0; u < n; ++u) {
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
    vector<int> input_vertexes(n);
    input_vertexes[0] = s;
    int num_input_vertexes = 1;;
	int current_p = p;
	//int current_p = min(p, num_input_vertexes);

    vector< vector<int> > prefix_sum_degrees(p);
    vector<int> prefix_sum_degrees_per_core;
    prefix_sum_degrees_per_core.resize(p);

    int * num_degrees = (int *) malloc(sizeof(int) * p);
    int * offset_degrees = (int *) malloc(sizeof(int) * p);
    int * num_vertexes = (int *) malloc(sizeof(int) * p);
    int * offset_vertexes = (int *) malloc(sizeof(int) * p);
    vector<int> queue_size(p);
    //int * queue_size = (int *) malloc(sizeof(int) * p);
    q.resize(p);
    while (num_input_vertexes > 0) {
		int new_p = min(p, num_input_vertexes);
		if(current_p != new_p){
			current_p = min(p, num_input_vertexes);
			prefix_sum_degrees_per_core.resize(current_p);
		}
		//current_p = p;
        _Cilk_for (int i = 0; i < current_p; i++) {
            get_even_split_size_and_offset(current_p, i, num_input_vertexes, num_vertexes[i], offset_vertexes[i]);

            // Generate prefix sum for out-degrees in core i's portion of input_vertexes
            prefix_sum_degrees[i].clear();
            prefix_sum_degrees[i].reserve(num_vertexes[i]);
            for (int j = 0; j < num_vertexes[i]; j++) {
                const int u = input_vertexes[j+offset_vertexes[i]];
                prefix_sum_degrees[i].push_back(adj[u].size());
            }
            destructive_serial_prefix_sum(prefix_sum_degrees[i]);
            assert((int)prefix_sum_degrees[i].size() == num_vertexes[i]);
            prefix_sum_degrees_per_core[i] = prefix_sum_degrees[i].empty() ? 0 : prefix_sum_degrees[i].back();
        }
        destructive_parallel_prefix_sum(prefix_sum_degrees_per_core);
        int total_degrees = prefix_sum_degrees_per_core.back();
		int A_start[current_p], A_end[current_p];
		int B[current_p];
		if(total_degrees > 0){
			A_start[0] = A_end[0] = -1;
			B[0] = -1;
			if(prefix_sum_degrees_per_core[0] > 0){
				A_start[0] = 0;
				A_end[0] = get_slot_number(current_p, total_degrees, prefix_sum_degrees_per_core[0]);
			}
			_Cilk_for(int i = 1; i < current_p; ++i){
				int prev_last = 0;
				int current_last = 0;
				A_start[i] = A_end[i] = B[i] = -1;
				if(prefix_sum_degrees_per_core[i] > prefix_sum_degrees_per_core[i-1]){
					prev_last = get_slot_number(current_p, total_degrees, prefix_sum_degrees_per_core[i-1]);
					current_last = get_slot_number(current_p, total_degrees, prefix_sum_degrees_per_core[i]);
					assert(prev_last <= current_p);
					assert(current_last <= current_p);
					assert(current_last >= prev_last);
					if(prev_last < current_last){
						A_start[i] = prev_last + 1;
						A_end[i] = current_last;
					}
				}
			}
			_Cilk_for(int i = 0; i < current_p; ++i){
				if(A_start[i] != -1){
					for(int j = A_start[i]; j <= A_end[i]; ++j){
						B[j] = i;
					}
				}
			}
		}
#if DEBUG_ARRAY
	cout<<"Testing Started"<<endl;
	for(int i = 0; i < current_p; ++i){
		cout<<i<<" "<<A_start[i]<<" "<<A_end[i]<<" "<<B[i]<<" "<<total_degrees<<endl;
		assert(B[i] != -1);
	}
	cout<<"Testing Ended"<<endl;
#endif

        _Cilk_for (int i = 0; i < current_p; i++) {
            get_even_split_size_and_offset(current_p, i, total_degrees, num_degrees[i], offset_degrees[i]);

            q[i].clear();
            if (num_degrees[i] > 0){
				assert(B[i] != -1);
                const int sublist = B[i];
                const int degrees_in_earlier_lists = (sublist==0) ? 0 : prefix_sum_degrees_per_core[sublist-1];
                const int vertex_in_sublist = find_index_in_prefix_sum(offset_degrees[i]+1 - degrees_in_earlier_lists, prefix_sum_degrees[sublist]);
                const int degrees_in_earlier_vertexes = (vertex_in_sublist==0) ? 0 : prefix_sum_degrees[sublist][vertex_in_sublist-1];

                const int starting_degree = offset_degrees[i] - degrees_in_earlier_lists - degrees_in_earlier_vertexes;
                assert(starting_degree >= 0);
                const int vertexes_in_earlier_sublists = offset_vertexes[sublist];
                const int starting_vertex = vertex_in_sublist + vertexes_in_earlier_sublists;
                assert(starting_degree < (int)adj[input_vertexes[starting_vertex]].size());

                int remaining = num_degrees[i];
                int vertex = starting_vertex;
                int degree = starting_degree;
                while (remaining > 0) {
                    int u = input_vertexes[vertex];
                    int limit = min(remaining + degree, (int)adj[u].size());
                    assert(limit > degree);
                    for (int j = degree; j < limit; j++) {
                        int v = adj[u][j];
                        if (d[v] == INFINITY) {
                            owner[v] = i;  // Dedup/duplicate expansion optimization
                            d[v] = d[u] + 1;
                            if (!adj[v].empty()) {
                                q[i].push_back(v);
                            }
                        }
                    }
                    remaining -= limit - degree;
                    degree = 0;
                    vertex++;
                }
                assert(remaining == 0);
            }
        }
        // Delete duplicate items from queues.
        _Cilk_for (int i = 0; i < current_p; i++) {
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
            //It can be moved up to prev _Cilk_for if no dedupping is done.
            queue_size[i] = q[i].size();
        }
        destructive_parallel_prefix_sum(queue_size, current_p);
        num_input_vertexes = queue_size[current_p-1];
        assert(num_input_vertexes <= n);
        _Cilk_for (int i = 0; i < current_p; i++) {
            const int offset = (i==0) ? 0 : queue_size[i-1];
            for (int j = 0; j < (int)q[i].size(); j++) {
                input_vertexes[offset+j] = q[i][j];
            }
        }
#if DEBUG_PRINT
	cout<<"Remaining "<<num_input_vertexes<<endl;
	for(int i = 0; i < current_p; ++i){
		cout<<q[i].size()<<endl;
		for(int j = 0; j < (int) q[i].size(); ++j){
			cout<<q[i][j]<<" ";
		}
		cout<<endl;
	}
#endif
    }
    free(num_degrees);
    free(offset_degrees);
    free(num_vertexes);
    free(offset_degrees);
    //free(queue_size);
    cilk::reducer_max<int> maxd(0);
    _Cilk_for (int u = 0; u < n; ++u) {
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
        fflush(stdout);
#if QUIT_EARLY 
        fprintf(stderr, "QUITTING\n");
        break;
#endif
    }
}
static void test_get_slot_number(void) {
    const int maxp = 10;
    int size[maxp];
    int offset[maxp];
    for (int p = 1; p < maxp; p++) {
        for (int total = 0; total < 100; total++) {
			if(p > total) continue;
            for (int i = 0; i < p; i++) {
                Graph::get_even_split_size_and_offset(p, i, total, size[i], offset[i]);
            }
			for(int i = 1; i <= total; ++i){
				int slot = Graph::get_slot_number(p, total, i);
				assert(offset[slot] <= i);
				assert(slot < p);
				if(slot < p-1){
					assert(i - 1  < offset[slot+1] );
				}
			}
        }
    }
}


// Verifies this splits everything evenly and properly.
static void test_get_even_split_size_and_offset(void) {
    const int maxp = 10;
    int size[maxp];
    int offset[maxp];
    for (int p = 1; p < maxp; p++) {
        for (int total = 0; total < 100; total++) {
            for (int i = 0; i < p; i++) {
                Graph::get_even_split_size_and_offset(p, i, total, size[i], offset[i]);
            }
            int size_first = size[0];
            int last_offset = offset[0];
            assert(size_first >= 0);
            assert(last_offset == 0);

            int last_size = size_first;
            int count = size_first;
            for (int i = 1; i < p; i++) {
                assert(size[i] >= 0);
                count += size[i];
                assert(last_size == size[i] || last_size - 1 == size[i]);
                last_size = size[i];
                assert(offset[i] == last_offset + size[i-1]);
                last_offset = offset[i];
            }
            assert(count == total);
        }
    }
}

// Verifies prefix sum works in serial and parallel
static void test_prefix_sum(void) {
    vector<uint64_t> sum(64);
    sum.resize(64, 0);
    sum[0] = 1;
    for (size_t i = 1; i < sum.size(); i++) {
        sum.at(i) = (sum.at(i-1) << 1) | 1;
    }
    vector<uint64_t> v_p(64);
    vector<uint64_t> v_s(64);
    for (int size = 0; size < 64; size++) {
        v_p.resize(size, 0);
        v_s.resize(size, 0);
        for (int i = 0; i < size; i++) {
            v_p.at(i) = v_s.at(i) = 1ULL<<i;
        }
        Graph::destructive_parallel_prefix_sum(v_p);
        Graph::destructive_serial_prefix_sum(v_s);
        for (int i = 0; i < size; i++) {
            assert(v_s.at(i) == sum.at(i));
            assert(v_p.at(i) == sum.at(i));
        }

    }
}

int cilk_main(int argc, char *argv[]) {
    assert(argc == 5);
    struct timespec start_t;
    struct timespec end_t;
    int error;
    error = clock_gettime(CLOCK_REALTIME, &start_t);
    assert(error==0);

    string file_name(argv[1]);
    bool opt_c = atoi(argv[2]);
    bool do_parallel = atoi(argv[3]);
    bool do_test = atoi(argv[4]);

    if (do_test) {
        fprintf(stderr, "TESTING INSTEAD!!\n");
        fflush(stderr);
        test_get_even_split_size_and_offset();
        test_prefix_sum();
		test_get_slot_number();
        exit(0);
    }
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

    fprintf(stderr, "Time taken: %s\n", argv[1]);
    fprintf(stderr, "Time taken: %f\n", end_d - start_d);
//    p.init(filename);
    return 0;
}

