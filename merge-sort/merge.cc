#include "util.h"

int main(void) {
    struct run_info<3> r(true, 0, 0, 0, nullptr);
        // mmap input
        // while input not empty
            // consider the first ~M items (actually find X)
            // sort in memory
            // output = tempfile()
            //  create output
            //  truncate to X
            //  mmap output
            // dump to temp file
            // store temp file name next_level_structure
    // while true
    //  delete current_structure
    //  current_structure = next_level
    //  create empty next_level_structure
        // start x-way merge
        //  if numfiles <= x
        //      output = outfile
        //  else
        //      output = tempfile()
        //  create output
        //  truncate to sum of all the files we're merging
        //  mmap output
        //  mmap first x temp files
        //  do x-way merge on the x temp files to the output
        //  munmap x temp files
        //  delete x temp files
        //  munmap output
        //  add output to next_level_structure
        //  if numfiles <= x return

        // Temp dir must be empty
        //  file names are actually element offsets
        //  file name can be (original order) element offset in %d_%d.temp
        // temp file name is of c bytes
        // Allocate 2 files of c * ceil(N/M) bytes
        //
        // max num files = N/M
}
