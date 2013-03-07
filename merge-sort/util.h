#pragma once


#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <stdint.h>
#include <fcntl.h>
#include <string.h>
#include <malloc.h>
#include <sys/mman.h>


//#include <sys/types.h>
#include <sys/stat.h>


//#include <sys/types.h>
//#include <sys/stat.h>

//TODO:
//  "Open file"
//      Options:
//          existing filename (input, output) XOR temp file level/number
//          bool create (create implies write, no create implies read)
//          uint64_t size_if_create (alternate uint64_t num_elements_if_create)
//          RETURN:
//              int error
//              struct run_info* {
//                  uint64_t size;
//                  uint64_t elements;
//                  int fd;
//                  auto_something<char*> filename;
//                  void* start; //mmap pointer
//                  bool is_temp;
//                  uint32_t run_level;
//                  uint64_t run_number;
//              }
//  "close file"
//          in: struct run_info*
//          bool delete_file <- assert delete_file true implies is_temp
//


template <uint32_t FIXED_SIZE = 100>
struct run_info {
    public:
    uint64_t size;
    uint64_t num_elements;
    int fd;
    char* filename;
    void *mmapped_address;
    bool is_temp;
    uint32_t run_level;
    uint64_t run_number;

    private:
    bool is_open;


    public:
    void close(bool is_delete) {
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

    ~run_info() {
        assert(!this->is_open);
        free(this->filename);
    }


    private:
    // open file (possibly create)
    run_info(bool is_create, bool is_temp, uint32_t run_level, uint64_t run_number, char* filename, uint64_t num_elements, char* temp_dir) :
        is_temp(is_temp),
        run_level(run_level),
        run_number(run_number)
    {
        int open_flags = O_NOATIME|O_DIRECT;
        if (is_create) {
            assert(num_elements > 0);
            this->num_elements = num_elements;
            this->size = this->num_elements * FIXED_SIZE;
            if (is_temp) {
                assert(!filename);
                // 8 hex for run_level, 16 hex for run_number, 1 for null
                if (temp_dir[strlen(temp_dir)-1] == '/') {
                    temp_dir[strlen(temp_dir)-1] = '\0';
                }
                ssize_t bufsize = strlen(temp_dir) + strlen("/") + strlen("_.temp") + 8 + 16 + 1;
                this->filename = (char*)malloc(bufsize);
                assert(this->filename);
                int c = snprintf(this->filename, bufsize, "%s/%08" PRIX32 "_%016" PRIX64 ".temp", temp_dir, run_level, run_number);
                assert(c == bufsize-1);
            } else {
                this->filename = strdup(filename);
                assert(this->filename);
            }
            open_flags |= O_WRONLY|O_CREAT|O_EXCL|O_TRUNC;
            this->fd = open(this->filename, open_flags, S_IWUSR|S_IRUSR);
            assert(this->fd >= 0);

            int r = ftruncate(this->fd, this->size);
            assert(r == 0);
        } else {
            assert(!num_elements);
            assert(filename);
            this->filename = strdup(filename);
            assert(this->filename);
            open_flags |= O_RDONLY;
            this->fd = open(this->filename, open_flags, S_IWUSR|S_IRUSR);
            assert(this->fd >= 0);
            struct stat buf;
            int r = fstat(this->fd, &buf);
            assert(r==0);
            this->size = buf.st_size;
            assert(this->size % FIXED_SIZE == 0);
            this->num_elements = this->size / FIXED_SIZE;
        }
        int prot_flags = is_create ? PROT_WRITE : PROT_READ;
        int option_flags = MAP_SHARED;
        // |  | MAP_HUGETLB;
        // | MAP_NORESERVE; (possible optimization)
        //
        this->size = this->num_elements * FIXED_SIZE;
        this->mmapped_address = mmap(nullptr, this->size, prot_flags, option_flags, this->fd, 0);
        assert(this->mmapped_address != MAP_FAILED);
        this->is_open = true;
    }

    public:
    run_info(bool is_create, uint32_t run_level, uint64_t run_number, uint64_t num_elements, char* temp_dir)
        : run_info(is_create, true, run_level, run_number, nullptr, num_elements, temp_dir) {}

    run_info(bool is_create, char* filename, uint64_t num_elements, char* temp_dir)
        : run_info(is_create, false, 0, 0, filename, num_elements, temp_dir) {}
};

