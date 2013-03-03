
#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include <cstdlib>
#include <iterator>
#include <vector>
#include <functional>
#include <limits>
#include <exception>
#include <iostream>

#include <funnel.h>
#include <sort.h>

using namespace std;
using namespace iosort;

class selfdescribing_exception : public exception
{
public:
	virtual std::ostream& describe(std::ostream& s) const = 0;
};

template<class T>
class bad_checksum_exception : public selfdescribing_exception
{
public:
	bad_checksum_exception(T was, T expected) : was(was), expected(expected) { }

	virtual const char *what() const throw()
		{ return "Bad checksum."; }
	virtual std::ostream& describe(std::ostream& s) const
		{ return s << "Expected checksum " << expected << ". Actual: " << was << '.' << endl; }
private:
	T was;
	T expected;
};

template<class T>
class out_of_order_exception : public selfdescribing_exception
{
public:
	out_of_order_exception(size_t one, size_t other, T first_value, T second_value) : one(one), other(other), first_value(first_value), second_value(second_value) { }

	virtual const char *what() const throw()
		{ return "Out of order."; }
	virtual std::ostream& describe(std::ostream& s) const
		{ return s << "Elements " << one << " and " << other << " out of order: " << first_value << '>' << second_value << endl; }
private:
	size_t one;
	size_t other;
	T first_value;
	T second_value;
};

template<class T, class Pred=less<T> >
class sorted_stream
{
	class iterator : public std::iterator<std::random_access_iterator_tag, T, int>
	{
		friend class sorted_stream<T>;
	public:
		inline iterator& operator++()
		{
			if( cnt )
			{
				if( s->pred(s->current,s->prev) )
					throw out_of_order_exception<T>(cnt, cnt+1, s->prev, s->current);
				s->prev = s->current;
				s->chksum += s->current;
			}
			else
			{
				s->prev = s->current;
				s->chksum = s->current;
			}
			++cnt;
			return *this;
		}
		inline iterator operator+(int i) const
		{
			iterator it = *this;
			it.cnt += i;
			return it;
		}
		inline T& operator*()
			{ return s->current; }
		inline const T& operator*() const
			{ return s->current; }
		inline bool operator==(const iterator& rhs) const
			{ return cnt == rhs.cnt; }
		inline bool operator!=(const iterator& rhs) const
			{ return cnt != rhs.cnt; }
		inline long operator-(const iterator& rhs) const
			{ return static_cast<long>(cnt-rhs.cnt); }
	public:
		iterator(sorted_stream *s) : s(s), cnt(0) { }
		iterator(sorted_stream *s, size_t c) : s(s), cnt(c) { }
	private:
		sorted_stream *s;
		size_t cnt;
	};

public:
	sorted_stream(size_t s) : size(s), pred() { }
	sorted_stream(size_t s, Pred pred) : size(s), pred(pred) { }
	iterator begin()
		{ return iterator(this); }
	iterator end()
		{ return iterator(this, size); }
	inline T checksum() const
		{ return chksum; }
protected:
	T prev;
	T current;
	T chksum;
	size_t size;
	Pred pred;
};


template<class value_type>
struct random_source
{
	value_type operator()();
	value_type checksum();
	void clear_sum();
};

template<>
struct random_source<int>
{
	struct generator_
	{
		generator_(random_source<int> *source) : source(source) { }
		int operator()()
		{
			int ret = rand();
			source->sum += ret;
			return ret;
		}
	private:
		random_source<int> *source;
	};

	generator_ generator()
	{
		return generator_(this);
	}
	random_source()
	{
		clear_sum();
	}
	int checksum()
	{
		return sum;
	}
	void clear_sum()
	{
		sum = 0;
	}
private:
	int sum;
};

class element
{
	friend std::ostream& operator<<(std::ostream& os, const element& e);
public:
	element() : key(), ptr() { };
	element(unsigned long int key, void* ptr) : key(key), ptr(ptr) { };
	inline bool operator<(const element& rhs) const
		{ return key<rhs.key; }
	inline bool operator==(const element& rhs) const
		{ return key==rhs.key; }
	inline bool operator!=(const element& rhs) const
		{ return key!=rhs.key; }

	inline element operator+=(const element& rhs)
		{ key+=rhs.key; return *this; }
private:
	unsigned long int key;
	void* ptr;
};
std::ostream& operator<<(std::ostream& os, const element& e)
{
	return os << '(' << e.key << ',' << e.ptr << ')';
}
namespace std
{
template<>
class numeric_limits<element>
{
public:
	inline static element min()
	{
		return element(std::numeric_limits<unsigned long int>::min(), static_cast<char*>(NULL)+0xF00DF00D);
	}
	inline static element max()
	{
		return element(std::numeric_limits<unsigned long int>::max(), static_cast<char*>(NULL)+0xF00DF00D);
	}
};
} // namespace std

template<>
struct random_source<element>
{
	struct generator_
	{
		generator_(random_source<element> *source) : source(source) { }
		element operator()()
		{
			element ret = element(rand(), reinterpret_cast<char*>(NULL)+0xBAADF00D);
			source->sum += ret;
			return ret;
		}
	private:
		random_source<element> *source;
	};

	generator_ generator()
	{
		return generator_(this);
	}

	random_source()
	{
		clear_sum();
	}
	element checksum()
	{
		return sum;
	}
	void clear_sum()
	{
		sum = element();
	}
private:
	element sum;
};


template<class Merger, class InpRanIt, class OutRanIt>
void check_merge_sort_(InpRanIt input_begin, InpRanIt input_end, OutRanIt output_begin, OutRanIt output_end)
{
	random_source<typename Merger::value_type> rand;
	generate(input_begin, input_end, rand.generator());

	cout << "Sorting " << input_end-input_begin << " elements..." << flush;
	merge_sort<Merger,typename Merger::splitter>(input_begin, input_end, output_begin);
	sorted_stream<typename Merger::value_type, typename Merger::predicate> ss(input_end-input_begin);
	copy(output_begin, output_end, ss.begin());
	if( ss.checksum() != rand.checksum() )
		throw bad_checksum_exception<typename Merger::value_type>(ss.checksum(), rand.checksum());

	cout << " and back again..." << flush;
	merge_sort<Merger,typename Merger::splitter>(output_begin, output_end, input_begin);
	sorted_stream<typename Merger::value_type, typename Merger::predicate> ssb(input_end-input_begin);
	copy(input_begin, input_end, ssb.begin());
	if( ssb.checksum() != rand.checksum() )
		throw bad_checksum_exception<typename Merger::value_type>(ss.checksum(), rand.checksum());

	cout << " done." << endl;
}

template<class Merger, class Container>
void check_merge_sort(size_t size)
{
	Container input(size);
	Container output(size);
	check_merge_sort_<Merger>(input.begin(), input.end(), output.begin(), output.end());
}

template<class Merger>
void check_merge_sort(size_t size)
{
	typename Merger::value_type* input = new typename Merger::value_type[size];
	typename Merger::value_type* output = new typename Merger::value_type[size];
	try
	{
		check_merge_sort_<Merger>(input, input+size, output, output+size);
	}
	catch( ... )
	{
		delete[] input;
		delete[] output;
		throw;
	}
	delete[] input;
	delete[] output;
}

void base_sort_test()
{
	// TODO: Implement
}

void funnel_test()
{
	// TODO: Implement
}

struct odd_before_even
{
	bool operator()(int lhs, int rhs)
	{
		if( (lhs&1) && !(rhs&1) )
			return true;
		if( !(lhs&1) && (rhs&1) )
			return false;
		return lhs < rhs;
	}
};

void merge_sort_test()
{
	for( size_t size=19; size!=19000; ++size )
	{
		check_merge_sort<merge_tree<int*,2> >(size);
		check_merge_sort<merge_tree<int*,4> >(size);
		check_merge_sort<merge_tree<element*,4> >(size);
		check_merge_sort<merge_tree<int*,2,default_splitter<2>,odd_before_even> >(size);
		check_merge_sort<merge_tree<int*,4,default_splitter<4>,odd_before_even> >(size);
		check_merge_sort<merge_tree<vector<int>::iterator,2>, vector<int> >(size);
		check_merge_sort<merge_tree<vector<int>::iterator,4>, vector<int> >(size);
		check_merge_sort<merge_tree<vector<element>::iterator,4>, vector<element> >(size);
		check_merge_sort<merge_tree<vector<int>::iterator,2,default_splitter<2>,odd_before_even>, vector<int> >(size);
		check_merge_sort<merge_tree<vector<int>::iterator,4,default_splitter<4>,odd_before_even>, vector<int> >(size);
	}
}

int main(int, char*[])
{
	try
	{
		base_sort_test();
		funnel_test();
		merge_sort_test();
	}
	catch( const selfdescribing_exception& ex )
	{
		ex.describe(cerr);
		return 1;
	}
	catch( const exception& ex )
	{
		cerr << "Exception cought: " << ex.what() << endl;
		return 1;
	}
	catch( const char *ex )
	{
		cerr << "Exception cought: " << ex << endl;
		return 1;
	}
	catch( ... )
	{
		cerr << "Exception cought." << endl;
		return 1;
	}
	return 0;
}
