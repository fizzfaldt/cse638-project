//Choose only one of option 1 or option 2
void pss(array, start, end) {
    if (end > start) {
        pssu(a, start, end);
        pssd(a, start, end, 0);
        sync; //Option 1
    }
}

int ppsu(array, start, end) {
    int size = end - start;
    if (size == 1) {
        return array[start];
    }
    int x = cilk_spawn ppsu(array, start, (start+end)/2);
    int y = ppsu(array, (start+end)/2, end);
    sync;
    a[end-1] = x+y;
    return = x+y;
}

void ppsd(array, start, end, ps) {
    int size = end - start;
    if (size == 1) {
        array[start] += ps;
        return;
    }
    int sum_left = array[(start+end)/2 -1];
    cilk_spawn ppsd(array, start, (start+end)/2, ps);
    ppsd(array, (start+end)/2, end, ps+sum_left);
    sync; //Option 2
}
