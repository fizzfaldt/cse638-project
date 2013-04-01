parallel_bfs(G,s,d) {
    n = G.V.size();


    parallel for u in [0...n) {
        d[u] = infinity;
    }
    d[s] = 0;
    p = get_num_cores();
    D = new Array[0..p)
    Q = new Array[0..p)
    S = new Array[0..p)
    I = new Array[0..1);
    while !I.isempty() {
        //Not worth even considering stealing
        // Generating D sublists and S arrays
        parallel for i in [0..p) {
            if i >= I.size mod p then
                    size = floor(I.size/p);
            else
                    size = ceil(I.size/p);
            startindex = i * floor(I.size/p) + min(i, I.size mod p);
            sum = 0;
            D[i] = new Array[0..size)
            for u in [startindex ... startindex + size) {
                sum += Gamma(I[u]).size;
                D[i][u-startindex] = sum;
            }
            S[i] = D[i][size-1]
        }
        //Not worth even considering stealing
        parallel_prefix_sum(S, p);
        W = S[p-1];
        //Maybe has room to add work stealing, can be imbalanced.
        //Running each work item (out degree)
        parallel for i in [0..p) {
            if i >= W mod p then
                    size = floor(W/p);
            else
                    size = ceil(W/p);
            startindex = i * floor(W/p) + min(i, W mod p);
            sublist = binary_search_index(startindex, S);
            vertex = binary_search_index(startindex - S[sublist-1], D[sublist]); // Treat S[-1] as 0
            degree = startindex - S[sublist-1] - D[sublist][vertex-1]; // Treat S[-1] and D[sublist][-1] as 0
            vertex += S[sublist-1]; //now an index into I[] //TODO: <_ this is blatantly wrong!! wrong array?

            Q[i] = new Queue
            {
                remaining = size
                while remaining > 0 {
                    u = I[vertex];
                    limit = min(remaining + degree, Gamma(u).size);
                    for v in Gamma(u)[degree...limit) {
                        if d[v] = infinity {
                            o[v] = i; // OPTIONAL (dedup/duplicate exploration)
                            d[v] = d[u] + 1;
                            q[i].enqueue(v);
                        }
                    }
                    remaining -= limit - degree;
                    degree = 0;
                    vertex++;
                }
            }
        }
        // OPTIONAL (dedup/duplicate exploration)
        parallel for i in [0..p) {
            //NOTE: This can be balanced better by doing a search
            temp = new Queue;
            while !q[i].empty {
                u = q[i].pop();
                if o[u] == i:
                    temp.enqueue(u)
            }
            delete q[i];
            q[i] = temp;
        }
        parallel for i in [0..p) {
            QS[i] = q[i].size
        }
        parallel_prefix_sum(QS, p);
        I = new Array[0..QS[p-1]);
        //Could be balanced better, alternatively make the for a pfor to use cilk stealing
        parallel for i in [0..p) {
            //NOTE: This can be balanced better by doing a search
            start = QS[i-1]; //treat QS[-1] as 0

            //NOTE: making this one parallel for can improve balance by using cilk-stealing
            for x in [0..q[i]size) {
                I[x+start] = q[x]
            }
        }
    }
}



