#pragma once

#define cilk_spawn
#define cilk_sync
#define cilk_for for
#define cilk_main main

//#define DISABLE_CILK
//#ifdef DISABLE_CILK
inline int __cilkrts_get_nworkers(void) {
    return 1;
}

namespace cilk {
    class context {
        public:
            int get_worker_count() { return 8; }
    };
    template<class T>
    class reducer_opand {
        private:
            T value;
        public:
            reducer_opand(bool initial) {
                value = initial;
            }

            reducer_opand & operator&=(const reducer_opand &rhs) {
                value &= rhs.value;
                return *this;
            }

            reducer_opand & operator&=(const T &rhs) {
                value &= rhs;
                return *this;
            }

            T get_value() {
                return value;
            }

    };

    template<class T>
    class reducer_max {
        private:
            T value;
        public:
            reducer_max(T value) {
                this->value = value;
            }

            reducer_max & operator=(const reducer_max &rhs) {
                value = rhs.value;
                return *this;
            }

            reducer_max & operator =(const T &rhs) {
                value = rhs;
                return *this;
            }

            T get_value() {
                return value;
            }
    };

    template<class T>
    T max_of(reducer_max<T> &x, T rhs) {
        return std::max(x.get_value(), rhs);
    }

    template<class T>
    class reducer_opadd {
        private:
            T value;
        public:
            reducer_opadd(void) {
                value = 0;
            }

            reducer_opadd & operator+=(const reducer_opadd &rhs) {
                value += rhs.value;
                return *this;
            }

            reducer_opadd & operator+=(const T &rhs) {
                value += rhs;
                return *this;
            }

            T get_value() {
                return value;
            }

    };
}

