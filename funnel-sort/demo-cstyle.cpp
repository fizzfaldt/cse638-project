
#include <funnel.h>
#include <sort.h>

#include <stdlib.h>

void sort_12k_cstyle()
{
	typedef int T;
	typedef iosort::merge_tree<T*,4> merger;

	size_t size = 12000;
	T *input = (T*)malloc(size*sizeof(T));
	T *output = (T*)malloc(size*sizeof(T));

	// populate input
	for( size_t i=0; i!=size; ++i )
		input[i] = rand();

	// sort
	iosort::merge_sort<merger,iosort::default_splitter<4> >(input, input+size, output);

	// output now contains copies of the elements in input in sorted order
}

int main(int argc, char* argv[])
{
	sort_12k_cstyle();
	return 0;
}
