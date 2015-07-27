#include <stdio.h>
#include <stdlib.h>
#include "htslib/vcf.h"
#include "htslib/khash.h"
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include "bcftools.h"

#ifdef DEBUG
#define ASSERT(X)  assert(X)
#else 
#define ASSERT(X)  ;
#endif

int int_comparator(const void* a, const void* b)
{
    int x = *((const int*)a);
    int y = *((const int*)b);
    return (x<y);
}

void destroy_PLmedian(plmedian_struct* info)
{
  if(info->m_output_fptr)
  {
    if(info->m_buffer_len)
      free(info->m_buffer);
    if(info->m_median_result_len)
      free(info->m_median_result);
    if(info->m_reorg_buffer_len)
      free(info->m_reorg_buffer);
    fclose(info->m_output_fptr);
  }
}

int compute_PLmedian(bcf_hdr_t* hdr, bcf1_t* line, int** median_result, int* median_result_len,
        int** buffer, int* buffer_len,
        int** reorg_buffer, int* reorg_buffer_len)
{
    bcf_unpack(line, BCF_UN_ALL);
    int num_vals = bcf_get_format_int32(hdr, line, "PL", buffer, buffer_len);
    if(num_vals < 0)
        return -1;
    int nsamples = bcf_hdr_nsamples(hdr);
    if(*reorg_buffer_len < nsamples)
    {
        *reorg_buffer = (int*)realloc((void*)*reorg_buffer, nsamples*sizeof(int));
        *reorg_buffer_len = nsamples;
    }
    int per_sample_veclen = num_vals/nsamples;
    if(*median_result_len < per_sample_veclen)
    {
        *median_result = (int*)realloc((void*)*median_result, per_sample_veclen*sizeof(int));
        *median_result_len = per_sample_veclen;
    }
    int i = 0;
#if 0
    for(i=0;i<nsamples;++i)
    {
        int j = 0;
        for(j=0;j<per_sample_veclen;++j)
            printf("%d,",(*buffer)[i*per_sample_veclen+j]);
        printf("\n");
    }
#endif
    for(i=0;i<per_sample_veclen;++i)
    {
        int j = 0;
        int offset = 0;
        int num_valid = 0;
        for(j=0;j<nsamples;++j)
        {
            int curr_value = (*buffer)[offset+i];
            if(curr_value != bcf_int32_missing && curr_value != bcf_int32_vector_end)
                (*reorg_buffer)[num_valid++] = curr_value;
            offset += per_sample_veclen;
        }
        if(num_valid > 0)
        {
            qsort(*reorg_buffer, num_valid, sizeof(int), int_comparator); 
            (*median_result)[i] = (*reorg_buffer)[num_valid/2];
        }
        else
            (*median_result)[i] = bcf_int32_missing;
    }
    return per_sample_veclen;
}

void print_PLmedian(FILE* fptr, bcf_hdr_t* hdr, bcf1_t* line, int* median_result, int median_result_len)
{
    fprintf(fptr,"%s,%d,%s,",bcf_seqname(hdr, line), line->pos+1, line->d.allele[0]);
    int i = 0;
    for(i=1;i<line->n_allele;++i)
    {
        if(i > 1)
            fprintf(fptr,",");
        fprintf(fptr, "%s", line->d.allele[i]);
    }
    fprintf(fptr, ",");
    for(i=0;i<median_result_len;++i)
    {
        if(i>0)
            fprintf(fptr,",");
        int curr_value = median_result[i];
        if(curr_value != bcf_int32_missing && curr_value != bcf_int32_vector_end)
            fprintf(fptr,"%d",curr_value);
        else
            fprintf(fptr,".");
    }
    fprintf(fptr,"\n");
}

int main_PLmedian(int argc, char *argv[])
{
    assert(argc >= 2 && "Requires at least 1 arguments <vcf1.gz>");
    //open file and read record
    htsFile* fptr = hts_open(argv[1], "r");
    assert(fptr);
    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    assert(hdr);
    bcf1_t* line = 0;
    int status = 0;
    line = bcf_init();
    
    int* median_result = 0;
    int median_result_len = 0;
    int* buffer = 0;
    int buffer_len = 0;
    int* reorg_buffer = 0;
    int reorg_buffer_len = 0;
    FILE* wfptr;
    if(argc >= 3)
      wfptr = fopen(argv[2],"w");
    else
      wfptr = stdout;
    int veclen = 0;
    while(1)
    {
        status = bcf_read(fptr, hdr, line);
        if(status < 0)
            break;
        veclen = compute_PLmedian(hdr, line, &median_result, &median_result_len,
                &buffer, &buffer_len,
                &reorg_buffer, &reorg_buffer_len);
        print_PLmedian(wfptr, hdr, line, median_result, veclen);
    }
    if(argc >= 3)
      fclose(wfptr);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fptr);
    return 1;
}

