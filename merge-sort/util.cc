#include "util.h"

template<uint32_t FIXED_SIZE>
void run_info<FIXED_SIZE>::init(bool is_delete) {

template<uint32_t FIXED_SIZE>
void run_info<FIXED_SIZE>::close(bool is_delete) {
    int r;
    assert(this->is_open);
    r = munmap(this->mmapped_address, this->size);
    assert(r==0);
    r = close(this->fd);
    assert(r==0);
    if (is_delete) {
        assert(this->is_temp);
        r = unlink(this->filename);
        assert(r==0);
    }
    this->is_open = false;
}

