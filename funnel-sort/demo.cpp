
#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include <functional>
#include <vector>
#include <cstdlib>

#include <funnel.h>
#include <sort.h>

using namespace std;
using namespace iosort;

void sort_12k()
{
	typedef merge_tree<vector<int>::iterator,4> merger;

	vector<int> input(12000);
	vector<int> output(12000);

	// populate input:
	generate(input.begin(), input.end(), rand);

	// sort
	merge_sort<merger,default_splitter<4> >(input.begin(), input.end(), output.begin());
 
	// output now contains copies of the elements in input in sorted order
}

int main(int, char*[])
{
	sort_12k();
	return 0;
}
