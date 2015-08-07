#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "bcftools.h"

main(int argc, char** argv)
{
    assert(argc >= 2 && "Needs <path_to_sqlite_file");
    void* ptr = allocate_sqlite3_mapping(argv[1]);
    //Given a row idx, get sample name
    //Space for sample_name must be allocated by caller
    char sample_name[100];
    query_sample_name(ptr, 0, sample_name);
    printf("Sample name for TileDB row idx 0 %s\n",sample_name);
    //Given contig names, return TileDB column idx for start of contig 
    const char* contig_names[] = { "1", "X" };
    const int64_t* contig_offset = query_contigs_offset(ptr, 2, contig_names, 0);
    printf("TileDB column ids for contigs 1 and X : %"PRIi64" and %"PRIi64"\n", contig_offset[0], contig_offset[1]);
    //Free internals of the structure
    free_sqlite3_data(ptr);
    //Free the structure
    free(ptr);
    return 0;
}
