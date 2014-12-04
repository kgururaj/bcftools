/*  vcfmerge.c -- Merge multiple VCF/BCF files to create one multi-sample file.

    Copyright (C) 2012-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <math.h>
#include <ctype.h>
#include "bcftools.h"
#include "vcmp.h"
/*For gVCFs reference is needed to get base pair at interval split locations*/
#include <htslib/faidx.h>

#include <htslib/khash.h>
KHASH_MAP_INIT_STR(strdict, int)
typedef khash_t(strdict) strdict_t;

#define SKIP_DONE 1
#define SKIP_DIFF 2

#define IS_VL_G(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_G)
#define IS_VL_A(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_A)
#define IS_VL_R(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_R)

#define USE_ID_MAP 1
// For merging INFO Number=A,G,R tags
typedef struct
{
    const char *hdr_tag;
    int type, nvals;
    int nbuf, mbuf;
    uint8_t *buf;
}
AGR_info_t;

// Rules for merging arbitrary INFO tags
typedef struct _info_rule_t
{
    char *hdr_tag;
    int merged_id;
    void (*merger)(bcf_hdr_t *hdr, bcf1_t *line, struct _info_rule_t *rule);
    int type;           // one of BCF_HT_*
    int block_size;     // number of values in a block
    int nblocks;        // number of blocks in nvals (the number of merged files)
    int nvals, mvals;   // used and total size of vals array
    void *vals;         // the info tag values
}
info_rule_t;

// Auxiliary merge data for selecting the right combination
//  of buffered records across multiple readers. maux1_t
//  corresponds to one buffered line.
typedef struct
{
    int skip;
    int *map;   // mapping from input alleles of a sample to the output array of the combined allele array
    int mmap;   // size of map array (only buffer[i].n_allele is actually used)
    int als_differ;
}
maux1_t;
typedef struct
{
    int n;  // number of readers
    char **als, **out_als;  // merged alleles (temp, may contain empty records) and merged alleles ready for output
    //nals : number of valid alleles in array als
    //mals : number of entries available in array als (some of them could be invalid)
    int nals, mals, nout_als, mout_als; // size of the output array
    int *cnt, ncnt; // number of records that refer to the alleles
    int *nbuf;      // readers have buffers of varying lengths
    int *smpl_ploidy, *smpl_nGsize; // ploidy and derived number of values in Number=G tags, updated for each line (todo: cache for missing cases)
    int *flt, mflt, minf;
    bcf_info_t *inf;// out_line's INFO fields
    bcf_fmt_t **fmt_map; // i-th output FORMAT field corresponds in j-th reader to i*nreader+j, first row is reserved for GT
    int nfmt_map;        // number of rows in the fmt_map array
    int *agr_map, nagr_map, magr_map;   // mapping between Number=AGR element indexes
    void *tmp_arr;
    int ntmp_arr;
    maux1_t **d;    // d[i][j] i-th reader, j-th buffer line
    AGR_info_t *AGR_info;
    int nAGR_info, mAGR_info;
    bcf_srs_t *files;
    int *has_line;  // which files are being merged
}
maux_t;

#ifdef COLLECT_STATS
enum StatEnum
{
    NUM_ACTIVE_SAMPLES_PER_LINE=0,
    NUM_FAST_SHAKE_BUFFER,
    LAST_STAT_IDX
};
char* StatNames[] = { "num_samples_per_line", "num_fast_shake_buffer", "NULL" };
typedef unsigned long long uint64;
typedef struct
{
    uint64* m_histogram;
    unsigned m_histogram_size;
    uint64 m_sum;
    double m_sum_square;
    uint64 m_count;
}stat_struct;
#endif

typedef struct
{
    int* m_id_2_merged_id;
    int m_num_ids;
    //Map from original to preprocessed
    int* m_original_id_2_preprocessed_id;
    //Map from original to preprocessed if the field needs to be copied somewhere else
    int* m_original_id_2_preprocessed_copy_id;
    int m_num_original_ids;
    //GQ id
    int m_GQ_id;
}reader_idmap;

typedef struct
{
    reader_idmap* m_readers_map;
    int** m_merged_id_2_reader_idmap; //2D map [nreaders][num_output_ids]
    int m_num_merged_ids;
    //map from id in merged hdr to index in merged bcf1_t* out
    int* m_merged_id_2_merged_idx; 
    //map from id in merged hdr to index in info_rules, -1 if none
    int* m_merged_id_2_info_rule_idx;
    //ids corresponding to DP
    int m_merged_DP_info_id;
    int m_merged_DP_format_id;
    //index vs/ id
    int m_merged_DP_info_idx;
    int m_merged_DP_format_idx;
    //array containing file ids with valid DP INFO fields in current record being merged
    int* m_DP_info_vals;
    //ids corresponding to MIN_DP
    int m_merged_MIN_DP_format_id;
    //index vs/ id
    int m_merged_MIN_DP_format_idx;
    //END id
    int m_merged_END_info_id;
    //flag whether merged line has only 1 variant of type NON_REF
    int m_merged_has_only_non_ref;
}idmap;

//calloc on args_t sets everything to 0
typedef struct
{
    FILE* m_output_fptr;
    int* m_median_result;
    int m_median_result_len;
    int* m_buffer;
    int m_buffer_len;
    int* m_reorg_buffer;
    int m_reorg_buffer_len;
}plmedian_struct;

int compute_PLmedian(bcf_hdr_t* hdr, bcf1_t* line, int** median_result, int* median_result_len,
        int** buffer, int* buffer_len,
        int** reorg_buffer, int* reorg_buffer_len);
void print_PLmedian(FILE* fptr, bcf_hdr_t* hdr, bcf1_t* line, int* median_result, int median_result_len);

typedef struct
{
    vcmp_t *vcmp;
    maux_t *maux;
    int header_only, collapse, output_type, force_samples, merge_by_id;
    char *header_fname, *output_fname, *regions_list, *info_rules, *file_list;
    info_rule_t *rules;
    int nrules;
    strdict_t *tmph;
    kstring_t tmps;
    bcf_srs_t *files;
    bcf1_t *out_line;
    htsFile *out_fh;
    bcf_hdr_t *out_hdr;
    //reference genome access - required for gVCF merge
    char* reference_filename;
    faidx_t* reference_faidx;
    char* reference_last_seq_read;
    int reference_last_read_pos;
    int reference_num_bases_read;
    char* reference_buffer;
    char* m_tag;
#ifdef COLLECT_STATS
    stat_struct stat_array[LAST_STAT_IDX];
#endif
    idmap m_idmap;
    char **argv;
    int argc;
    plmedian_struct  m_plmedian_info;
}
args_t;

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

kstring_t g_debug_string = { 0, 0, 0 };
#ifdef DEBUG
#define ASSERT(X)  assert(X)
FILE* g_debug_fptr = 0;
FILE* g_vcf_debug_fptr = 0;
#else
#define ASSERT(X) ;
#endif
unsigned g_preprocess_vcfs = 0;
unsigned g_is_input_gvcf = 0;
unsigned g_do_gatk_merge = 0;
const char* g_info2format_suffix = "_INFO";
unsigned g_measure_iterator_timing_only = 0;

//Merge configuration
enum MergeActionEnum
{
  MERGE_DEFAULT=0,
  MERGE_INFO,
  MERGE_FORMAT,
  MERGE_KEEP,
  MERGE_DROP,
  MERGE_MOVE_TO_FORMAT,
  MERGE_COPY_TO_FORMAT,
};
#define do_detach(action_type) ((action_type) == MERGE_DROP || (action_type) == MERGE_MOVE_TO_FORMAT)
typedef struct
{
    int is_initialized;
    int default_action;   //one of the values in MergeActionEnum
    strdict_t* info_field2action;   //map/dict from info field name to action
    strdict_t* format_field2action;   //map/dict from info field name to action
    int** file2id2action;        //2D array, 1 array for each reader, maps id in BCF_DT_ID dictionary to action
    int num_files;
}merge_config;
strdict_t*  g_merge_config_token_2_idx = 0;
merge_config g_merge_config = { 0, 0, 0, 0, 0, 0 };

int get_merge_action(const char* key_string, int is_info, int is_format)
{
    if(g_merge_config.is_initialized == 0)
        return MERGE_KEEP;
    khiter_t kitr = is_info ? kh_get(strdict, g_merge_config.info_field2action, key_string)
        : kh_get(strdict, g_merge_config.format_field2action, key_string);
    int action_type =
        ((is_info && kitr == kh_end(g_merge_config.info_field2action)) 
            || (is_format &&  kitr == kh_end(g_merge_config.format_field2action)))
        ? g_merge_config.default_action //if no entry found, do default action
        : is_info ? kh_val(g_merge_config.info_field2action, kitr) : kh_val(g_merge_config.format_field2action, kitr);
    return action_type;
}

#define REALLOC_IF_NEEDED(ptr, curr_num_elements, new_num_elements, type)       \
    if((new_num_elements) > (curr_num_elements))                                \
    {                                                                           \
        (curr_num_elements) = (new_num_elements);                               \
        (ptr) = (type*)realloc((ptr), (new_num_elements)*sizeof(type));         \
    }

#ifdef COLLECT_STATS
void initialize_stats(args_t* args)
{
    memset(args->stat_array, 0, LAST_STAT_IDX*sizeof(stat_struct));
    args->stat_array[NUM_ACTIVE_SAMPLES_PER_LINE].m_histogram = 
        (uint64*)calloc(args->files->nreaders+1, sizeof(uint64)); //0 to nreaders
    args->stat_array[NUM_ACTIVE_SAMPLES_PER_LINE].m_histogram_size = args->files->nreaders+1;
}

void destroy_stats(args_t* args)
{
    int i = 0;
    for(i=0;i<LAST_STAT_IDX;++i)
        if(args->stat_array[i].m_histogram)
            free(args->stat_array[i].m_histogram);
}

void update_stat(args_t* args, int stat_idx, int value)
{
    ASSERT(stat_idx >= 0 && stat_idx < LAST_STAT_IDX);
    stat_struct* ptr = &(args->stat_array[stat_idx]);
    ptr->m_sum += ((uint64)value);
    double v = (double)value;
    ptr->m_sum_square += v*v;
    ++(ptr->m_count);
    switch(stat_idx)
    {
        case NUM_ACTIVE_SAMPLES_PER_LINE:
            ASSERT(((unsigned)value) < ptr->m_histogram_size);
            ++(ptr->m_histogram[value]);
            break;
        default:
            break;
    }
}

void print_stats(args_t* args)
{
    char filename[100];
    if(args->m_tag)
        sprintf(filename,"%s_stats.csv",args->m_tag);
    else
        sprintf(filename,"stats.csv");
    FILE* fptr = fopen(filename,"w");
    int i = 0;
    for(i=0;i<LAST_STAT_IDX;++i)
    {
        stat_struct* ptr = &(args->stat_array[i]);
        ptr->m_count = (ptr->m_count == 0) ? 1 : ptr->m_count;
        double mean =  ((double)(ptr->m_sum))/(ptr->m_count);
        double variance = (ptr->m_sum_square/ptr->m_count) - (mean*mean);
        fprintf(fptr,"%s,%llu,%llu,%lf,%lf,%lf\n",StatNames[i], ptr->m_count, ptr->m_sum, mean,
                ptr->m_sum_square, variance);
        if(ptr->m_histogram)
        {
            int j = 0;
            for(j=0;j<ptr->m_histogram_size;++j)
                fprintf(fptr,"%d,%llu\n",j,ptr->m_histogram[j]);
        }
    }
    fclose(fptr);
}
#endif

/**
 * Unfortunately, this has to be a macro, because different data types are used, C++ templates would solve this
 * line: current reader's line being considered for merging
 * length_type: one of BCF_VL_*
 * als: maux1_t* which contains mapping for alleles
 * nout_als: number of alleles (including REF) in the merged line
 * src: correctly typecast array for the current field for the current reader's line
 * dst: correctly typecast array for the current field for the merged line
 * dst_type_t: type of the destination array
**/
#define initialize_ARG_vectors_for_NON_REF(line, length_type, als, nout_als, src, src_length, dst, dst_type_t) \
{ \
    /*FIXME: Works only for haploid and diploid*/ \
    /*copy values corresponding to NON_REF allele instead of missing value*/ \
    ASSERT(line->m_non_ref_idx >= 1 && line->m_non_ref_idx < line->n_allele); \
    ASSERT(line->d.var[line->m_non_ref_idx].type == VCF_NON_REF); \
    int l = 0; \
    if(length_type == BCF_VL_A || length_type == BCF_VL_R) \
    { \
        int exclude_ref = (length_type == BCF_VL_A) ? 1 : 0; \
        /*-exclude_ref because m_non_ref_idx includes REF also, however, vector of length A does not*/ \
        int ndst = nout_als - exclude_ref; \
        /*printf("Out num alleles %d input num allele %d non_ref idx %d\n", nout_als, line->n_allele, line->m_non_ref_idx);*/ \
        /*-exclude_ref because m_non_ref_idx includes REF also, however, vector of length A does not*/ \
        dst_type_t non_ref_val = (dst_type_t)(src[line->m_non_ref_idx-exclude_ref]); \
        /*printf("Non ref val %d\n",(int)non_ref_val);*/ \
        for(l=0;l<ndst;++l) \
            dst[l] = non_ref_val; \
    } \
    else \
    { \
        if(length_type == BCF_VL_G) \
        { \
            int is_diploid = (src_length > line->n_allele); \
            if(is_diploid) \
            { \
                /*ref|ref GT val*/ \
                dst[bcf_alleles2gt(0, 0)] = (dst_type_t)(src[bcf_alleles2gt(0,0)]); \
                /*values for ref|non_ref from src - copy for every GT where ref|<other> occurs*/ \
                dst_type_t ref_non_ref_gt_val = (dst_type_t)(src[bcf_alleles2gt(0, line->m_non_ref_idx)]); \
                /*values for non_ref|non_ref from src - copy for every GT without ref*/ \
                dst_type_t non_ref_non_ref_gt_val = (dst_type_t)(src[bcf_alleles2gt(line->m_non_ref_idx, line->m_non_ref_idx)]); \
                /*printf("GT format ref/non-ref %d non-ref/non-ref %d\n",(int)ref_non_ref_gt_val, (int)non_ref_non_ref_gt_val);*/ \
                for(l=1;l<nout_als;++l) \
                { \
                    /*Copy ref|non-ref GT values*/ \
                    dst[bcf_alleles2gt(0, l)] = ref_non_ref_gt_val; \
                    int m = 1; \
                    /*Copy non-ref|non-ref GT values*/ \
                    for(m=1;m<=l;++m) \
                        dst[bcf_alleles2gt(m, l)] = non_ref_non_ref_gt_val; \
                } \
                /*for GT of type var|<other> where var exists in current file, but other may or may not exist*/ \
                /*fill such entries with value corresponding to GT var|<non_ref>*/ \
                for(l=1;l<line->n_allele;++l) \
                { \
                    int lnew = als->map[l]; \
                    dst_type_t val = (dst_type_t)(src[bcf_alleles2gt(l, line->m_non_ref_idx)]); \
                    int m = 1; \
                    for(m=1;m<nout_als;++m) \
                        dst[bcf_alleles2gt(m, lnew)] = val; \
                } \
            } \
            else /*haploid*/ \
            { \
                dst_type_t non_ref_val = (dst_type_t)(src[line->m_non_ref_idx]); \
                for (l=0; l<nout_als; l++) { dst[l] = non_ref_val; } \
            } \
        } \
    } \
}

//Templates!
#define comparator(data_type) \
int data_type ## _comparator(const void* aptr, const void* bptr) \
{ \
    data_type a = *((data_type*)aptr); \
    data_type b = *((data_type*)bptr); \
    return (a < b) ? -1 : (a == b) ? 0 : 1; \
}

comparator(int32_t) 
comparator(float) 

static void info_rules_merge_median(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, ndim = rule->block_size;
    assert(ndim == 1 && rule->nvals == rule->nblocks && "Only single value fields are handled by median merge");
    #define BRANCH(type_t,is_missing, comparator) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        qsort(ptr, rule->nvals, sizeof(type_t), comparator); /*FIXME: use selection algorithm for median*/ \
        ptr[0] = ptr[rule->nvals/2]; /*middle element*/ \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, int32_t_comparator); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), float_comparator); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_sum(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) ptr[j] += ptr[j+i*ndim]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_avg(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (j=0; j<ndim; j++) \
        { \
            double sum = 0; \
            for (i=0; i<rule->nblocks; i++) sum += ptr[j+i*ndim]; \
            ptr[j] = sum / rule->nblocks; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_min(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] > ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MAX); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_max(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] < ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MIN); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), -HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_join(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    if ( rule->type==BCF_HT_STR )
    {
        ((char*)rule->vals)[rule->nvals] = 0;
        bcf_update_info_string(hdr,line,rule->hdr_tag,rule->vals);
    }
    else
        bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,rule->nvals,rule->type);
}

#ifndef USE_ID_MAP
static int info_rules_comp_key2(const void *a, const void *b)
{
    info_rule_t *rule1 = (info_rule_t*) a;
    info_rule_t *rule2 = (info_rule_t*) b;
    return strcmp(rule1->hdr_tag, rule2->hdr_tag);
}

static int info_rules_comp_key(const void *a, const void *b)
{
    char *key = (char*) a;
    info_rule_t *rule = (info_rule_t*) b;
    return strcmp(key, rule->hdr_tag);
}
#endif

static void info_rules_init(args_t *args)
{
    if ( args->info_rules && !strcmp("-",args->info_rules) ) return;

    kstring_t str = {0,0,0};
    if ( !args->info_rules )
    {
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP")) ) kputs("DP:sum",&str);
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP4")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("DP4:sum",&str);
        }
        if ( !str.l ) return;
        args->info_rules = str.s;
    }

    args->nrules = 1;
    char *ss = strdup(args->info_rules), *tmp = ss;
    int n = 0;
    while ( *ss )
    {
        if ( *ss==':' ) { *ss = 0; n++; if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        else if ( *ss==',' ) { *ss = 0; args->nrules++; n++; if ( n%2==1 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        ss++;
    }
    if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules);
    args->rules = (info_rule_t*) calloc(args->nrules,sizeof(info_rule_t));

    n = 0;
    ss = tmp;
    while ( n < args->nrules )
    {
        info_rule_t *rule = &args->rules[n];
        rule->hdr_tag = strdup(ss);
        int id = bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, rule->hdr_tag);
	if(id < 0)
	{
	  fprintf(stderr,"WARNING: output header does not contain tag %s - dropping associated merge rule\n",rule->hdr_tag);
	  ++n;
	  continue;
	}
	ASSERT(id >=0 && id < args->out_hdr->n[BCF_DT_ID]);
        ASSERT(strcmp(args->out_hdr->id[BCF_DT_ID][id].key, rule->hdr_tag) == 0);
        if ( !bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,id) ) error("The tag is not defined in the header: \"%s\"\n", rule->hdr_tag);
        rule->type = bcf_hdr_id2type(args->out_hdr,BCF_HL_INFO,id);
        if ( rule->type!=BCF_HT_INT && rule->type!=BCF_HT_REAL && rule->type!=BCF_HT_STR ) error("The type is not supported: \"%s\"\n", rule->hdr_tag);
#ifdef USE_ID_MAP
        ASSERT(id >= 0 && id < args->out_hdr->n[BCF_DT_ID]);
        args->m_idmap.m_merged_id_2_info_rule_idx[id] = n;
        rule->merged_id = id;
#endif

        while ( *ss ) ss++; ss++;
        if ( !*ss ) error("Could not parse INFO rules, missing logic of \"%s\"\n", rule->hdr_tag);

        int is_join = 0;
        if ( !strcasecmp(ss,"sum") ) rule->merger = info_rules_merge_sum;
        else if ( !strcasecmp(ss,"avg") ) rule->merger = info_rules_merge_avg;
        else if ( !strcasecmp(ss,"median") ) rule->merger = info_rules_merge_median;
        else if ( !strcasecmp(ss,"min") ) rule->merger = info_rules_merge_min;
        else if ( !strcasecmp(ss,"max") ) rule->merger = info_rules_merge_max;
        else if ( !strcasecmp(ss,"join") ) { rule->merger = info_rules_merge_join; is_join = 1; }
        else error("The rule logic \"%s\" not recognised\n", ss);

        if ( !is_join && rule->type==BCF_HT_STR )
            error("Numeric operation \"%s\" requested on non-numeric field: %s\n", ss, rule->hdr_tag);
        if ( bcf_hdr_id2number(args->out_hdr,BCF_HL_INFO,id)==0xfffff )
        {
            int is_agr = (
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_A ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_G ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_R
                    ) ? 1 : 0;
            if ( is_join && is_agr )
                error("Cannot -i %s:join on Number=[AGR] tags is not supported.\n", rule->hdr_tag);
            if ( !is_join && !is_agr )
                error("Only fixed-length vectors are supported with -i %s:%s\n", ss, rule->hdr_tag);
        }

        while ( *ss ) ss++; ss++; n++;
    }
    free(str.s);
    free(tmp);
#ifndef USE_ID_MAP
    qsort(args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key2);
#endif
}
static void info_rules_destroy(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
    {
        info_rule_t *rule = &args->rules[i];
        free(rule->hdr_tag);
        free(rule->vals);
    }
    free(args->rules);
}
static void info_rules_reset(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
        args->rules[i].nblocks = args->rules[i].nvals = args->rules[i].block_size = 0;
}
static int info_rules_add_values(args_t *args, bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule, maux1_t *als, int var_len)
{
    int ret = bcf_get_info_values(hdr, line, rule->hdr_tag, &args->maux->tmp_arr, &args->maux->ntmp_arr, rule->type);
    if ( ret<=0 )
        error("FIXME: error parsing %s at %s:%d .. %d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1,ret);
    rule->nblocks++;

    if ( rule->type==BCF_HT_STR )
    {
        int need_comma = rule->nblocks==1 ? 0 : 1;
        hts_expand(char,rule->nvals+ret+need_comma+1,rule->mvals,rule->vals);            // 1 for null-termination
        char *tmp = (char*) rule->vals + rule->nvals;
        if ( rule->nvals>0 ) { *tmp = ','; tmp++; }
        strncpy(tmp,(char*)args->maux->tmp_arr,ret);
        rule->nvals += ret + need_comma;
        return 1;
    }

    int i, j;
    if ( var_len==BCF_VL_A )
    {
        assert( ret==line->n_allele-1 );
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        // create mapping from source file ALT indexes to dst file indexes
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i+1] - 1;
        rule->block_size = args->maux->nout_als - 1;
    }
    else if ( var_len==BCF_VL_R )
    {
        assert( ret==line->n_allele );
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
        rule->block_size = args->maux->nout_als;
    }
    else if ( var_len==BCF_VL_G )
    {
        args->maux->nagr_map = bcf_alleles2gt(line->n_allele-1,line->n_allele-1)+1;
        assert( ret==line->n_allele || ret==args->maux->nagr_map );
        if ( ret==line->n_allele ) // haploid
        {
            args->maux->nagr_map = line->n_allele;
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
            rule->block_size = args->maux->nout_als;
        }
        else
        {
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            int k_src = 0;
            for (i=0; i<line->n_allele; i++)
            {
                for (j=0; j<=i; j++)
                {
                    args->maux->agr_map[k_src] = bcf_alleles2gt(als->map[i],als->map[j]);
                    k_src++;
                }
            }
            rule->block_size = bcf_alleles2gt(args->maux->nout_als-1,args->maux->nout_als-1)+1;
        }
    }
    else
    {
        //FIXME: will not work for non-deterministic variable length INFO fields
        if ( rule->nblocks>1 && ret!=rule->block_size )
            error("Mismatch in number of values for INFO/%s at %s:%d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1);
        rule->block_size = ret;
        args->maux->nagr_map = 0;
    }
    maux_t* ma = args->maux;

    #define BRANCH(src_type_t,dst_type_t,set_missing) { \
        src_type_t *src = (src_type_t *) args->maux->tmp_arr; \
        hts_expand0(dst_type_t,(rule->nvals+rule->block_size),rule->mvals,rule->vals); \
        dst_type_t *dst = (dst_type_t *) rule->vals + rule->nvals; \
        rule->nvals += rule->block_size; \
        if ( !args->maux->nagr_map ) \
        { \
            for (i=0; i<ret; i++) dst[i] = src[i]; \
        } \
        else \
        { \
            if(line->m_non_ref_idx >= 0) \
            { \
                initialize_ARG_vectors_for_NON_REF(line, var_len, als, (ma->nout_als), src, ret, dst, dst_type_t); \
            } \
            else \
                for (i=0; i<rule->block_size; i++) set_missing; \
            for (i=0; i<ret; i++) dst[args->maux->agr_map[i]] = src[i]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int, int32_t, dst[i] = bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, float, bcf_float_set_missing(dst[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    return 1;
}

int bcf_hdr_sync(bcf_hdr_t *h);

void bcf_hdr_merge(bcf_hdr_t *hw, const bcf_hdr_t *hr, const char *clash_prefix, int force_samples)
{
    // header lines
    int ret = bcf_hdr_combine(hw, hr);
    if ( ret!=0 ) error("Error occurred while merging the headers.\n");

    // samples
    int i;
    for (i=0; i<bcf_hdr_nsamples(hr); i++)
    {
        char *name = hr->samples[i];
        if ( bcf_hdr_id2int(hw, BCF_DT_SAMPLE, name)!=-1 )
        {
            // there is a sample with the same name
            if ( !force_samples ) error("Error: Duplicate sample names (%s), use --force-samples to proceed anyway.\n", name);

            int len = strlen(hr->samples[i]) + strlen(clash_prefix) + 1;
            name = (char*) malloc(sizeof(char)*(len+1));
            sprintf(name,"%s:%s",clash_prefix,hr->samples[i]);
            bcf_hdr_add_sample(hw,name);
            free(name);
        }
        else
            bcf_hdr_add_sample(hw,name);
    }
}

void debug_als(char **als, int nals)
{
    int k; for (k=0; k<nals; k++) fprintf(stderr,"%s ", als[k]);
    fprintf(stderr,"\n");
}

/**
 * normalize_alleles() - create smallest possible representation of the alleles
 * @als:    alleles to be merged, first is REF (rw)
 * @nals:   number of $a alleles
 *
 * Best explained on an example:
 *      In:  REF=GTTT  ALT=GTT
 *      Out: REF=GT    ALT=G
 *
 * Note: the als array will be modified
 */
void normalize_alleles(char **als, int nals)
{
    if ( !als[0][1] ) return;   // ref is 1base long, we're done

    int j, i = 1, done = 0;
    int *lens = (int*) malloc(sizeof(int)*nals);
    for (j=0; j<nals; j++) lens[j] = strlen(als[j]);

    while ( i<lens[0] )
    {
        for (j=1; j<nals; j++)
        {
            if ( i>=lens[j] ) done = 1;
            if ( als[j][lens[j]-i] != als[0][lens[0]-i] ) { done = 1; break; }
        }
        if ( done ) break;
        i++;
    }
    if ( i>1 )
    {
        i--;
        als[0][lens[0]-i] = 0;
        for (j=1; j<nals; j++) als[j][lens[j]-i] = 0;
    }
    free(lens);
}

 /**
 * merge_alleles() - merge two REF,ALT records, $a and $b into $b.
 * @a:      alleles to be merged, first is REF
 * @na:     number of $a alleles
 * @map:    map from the original $a indexes to new $b indexes (0-based)
 * @b:      alleles to be merged, the array will be expanded as required
 * @nb:     number of $b alleles (ptr because this int will be changed)
 * @mb:     size of $b
 * @is_split_record: A terrible hack, flag set to 1, if the BCF record for this allele was produced
 *          by a gVCF interval split (I'll stand in the corner for now)
 *
 * Returns NULL on error or $b expanded to incorporate $a alleles and sets
 * $map. Best explained on an example:
 *      In:     REF   ALT
 *           a: ACG,  AC,A    (1bp and 2bp deletion)
 *           b: ACGT, A       (3bp deletion)
 *      Out:
 *           b: ACGT, A,ACT,AT (3bp, 1bp and 2bp deletion)
 *           map: 0,2,3
 * Here the mapping from the original $a alleles to the new $b alleles is 0->0,
 * 1->2, and 2->3.
 */
char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb, variant_t* var_types, int* contains_non_ref, int is_split_record)
{
    // reference allele never changes
    map[0] = 0;

    int i,j;
    int rla = !a[0][1] ? 1 : strlen(a[0]);
    int rlb = !b[0][1] ? 1 : strlen(b[0]);

#ifdef USE_SPLIT_RECORD_HACK
    if(is_split_record)
        rla = 1;
#endif

    // the most common case: same SNPs
    if ( na==2 && *nb==2 && rla==1 && rlb==1 && a[1][0]==b[1][0] && !a[1][1] && !b[1][1] )
    {
        map[1] = 1;
        return b;
    }

#ifdef USE_SPLIT_RECORD_HACK
    if(is_split_record)
    {
        a[0][0] = b[0][0];
        a[0][1] = '\0';
    }
#endif

    // Sanity check: reference prefixes must be identical
    if ( strncmp(a[0],b[0],rla<rlb?rla:rlb) )
    {
        if ( strncasecmp(a[0],b[0],rla<rlb?rla:rlb) )
        {
            fprintf(stderr, "The REF prefixes differ: %s vs %s (%d,%d)\n", a[0],b[0],rla,rlb);
            return NULL;
        }
        // Different case, change to uppercase
        for (i=0; i<na; i++)
        {
            int len = strlen(a[i]);
            for (j=0; j<len; j++) a[i][j] = toupper(a[i][j]);
        }
        for (i=0; i<*nb; i++)
        {
            int len = strlen(b[i]);
            for (j=0; j<len; j++) b[i][j] = toupper(b[i][j]);
        }
    }

    int n = *nb + na;
    hts_expand0(char*,n,*mb,b);

    // $b alleles need expanding
    if ( rla>rlb )
    {
        for (i=0; i<*nb; i++)
        {
            int l = strlen(b[i]);
            b[i] = (char*) realloc(b[i],l+rla-rlb+1);
            memcpy(b[i]+l,a[0]+rlb,rla-rlb+1);
        }
    }

    // now check if the $a alleles are present and if not add them
    for (i=1; i<na; i++)
    {
        char *ai;
        if ( rlb>rla )  // $a alleles need expanding
        {
            int l = strlen(a[i]);
            ai = (char*) malloc(l+rlb-rla+1);
            memcpy(ai,a[i],l);
            memcpy(ai+l,b[0]+rla,rlb-rla+1);
        }
        else
            ai = a[i];

        //handle NON_REFs later to ensure it's always the last of the list
        if(var_types[i].type & VCF_NON_REF)
        {
            *contains_non_ref = 1;
            if ( rlb>rla ) free(ai);
            continue;
        }

        for (j=1; j<*nb; j++)
            if ( !strcmp(ai,b[j]) ) break;

        if ( j<*nb ) // $b already has the same allele
        {
            map[i] = j;
            if ( rlb>rla ) free(ai);
            continue;
        }
        // new allele
        map[i] = *nb;
        b[*nb] = rlb>rla ? ai : strdup(ai);
        (*nb)++;
    }
    return b;
}

maux_t *maux_init(bcf_srs_t *files)
{
    maux_t *ma = (maux_t*) calloc(1,sizeof(maux_t));
    ma->n      = files->nreaders;
    ma->nbuf   = (int *) calloc(ma->n,sizeof(int));
    ma->d      = (maux1_t**) calloc(ma->n,sizeof(maux1_t*));
    ma->files  = files;
    int i, n_smpl = 0;
    for (i=0; i<ma->n; i++)
        n_smpl += bcf_hdr_nsamples(files->readers[i].header);
    ma->smpl_ploidy = (int*) calloc(n_smpl,sizeof(int));
    ma->smpl_nGsize = (int*) malloc(n_smpl*sizeof(int));
    ma->has_line = (int*) malloc(ma->n*sizeof(int));
    return ma;
}
void maux_destroy(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) // for each reader
    {
        if ( !ma->d[i] ) continue;
        int j;
        for (j=0; j<ma->nbuf[i]; j++)  // for each buffered line
            if ( ma->d[i][j].map ) free(ma->d[i][j].map);
        free(ma->d[i]);
    }
    for (i=0; i<ma->mAGR_info; i++)
        free(ma->AGR_info[i].buf);
    free(ma->agr_map);
    free(ma->AGR_info);
    if (ma->ntmp_arr) free(ma->tmp_arr);
    if (ma->nfmt_map) free(ma->fmt_map);
    // ma->inf freed in bcf_destroy1
    free(ma->d);
    free(ma->nbuf);
    for (i=0; i<ma->mals; i++) free(ma->als[i]);
    if (ma->mout_als) free(ma->out_als);
    free(ma->als);
    free(ma->cnt);
    free(ma->smpl_ploidy);
    free(ma->smpl_nGsize);
    free(ma->has_line);
    free(ma);
}
void maux_expand1(maux_t *ma, int i)
{
    if ( ma->nbuf[i] <= ma->files->readers[i].last_valid_record_idx )
    {
        int n = ma->files->readers[i].last_valid_record_idx + 1;
        ma->d[i] = (maux1_t*) realloc(ma->d[i], sizeof(maux1_t)*n);
        memset(ma->d[i]+ma->nbuf[i],0,sizeof(maux1_t)*(n-ma->nbuf[i]));
        ma->nbuf[i] = n;
    }
}
void maux_reset(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) maux_expand1(ma, i);
    for (i=1; i<ma->ncnt; i++) ma->cnt[i] = 0;
}
void maux_debug(maux_t *ma, int ir, int ib)
{
    printf("[%d,%d]\t", ir,ib);
    int i;
    for (i=0; i<ma->nals; i++)
    {
        printf(" %s [%d]", ma->als[i], ma->cnt[i]);
    }
    printf("\n");
}

void merge_chrom2qual(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    kstring_t *tmps = &args->tmps;
    tmps->l = 0;

    maux_t *ma = args->maux;
    int *al_idxs = (int*) calloc(ma->nals,sizeof(int));
    bcf_float_set_missing(out->qual);

    // CHROM, POS, ID, QUAL
    out->pos = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;       //set in merge_buffer

        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->processed_header;

        // alleles
        int j;
        for (j=1; j<line->n_allele; j++)
            al_idxs[ ma->d[i][0].map[j] ] = 1;

        // position
        if ( out->pos==-1 )
        {
            const char *chr = hdr->id[BCF_DT_CTG][line->rid].key;
            out->rid = bcf_hdr_name2id(out_hdr, chr);
            if ( strcmp(chr,out_hdr->id[BCF_DT_CTG][out->rid].key) ) error("Uh\n");
            out->pos = line->pos;
        }

        // ID
        if ( line->d.id[0]!='.' || line->d.id[1] )
        {
            kitr = kh_get(strdict, tmph, line->d.id);
            if ( kitr == kh_end(tmph) )
            {
                if ( tmps->l ) kputc(';', tmps);
                kputs(line->d.id, tmps);
                kh_put(strdict, tmph, line->d.id, &ret);
            }
        }
        // set QUAL to the max qual value. Not exactly correct, but good enough for now
        if ( !bcf_float_is_missing(files->readers[i].buffer[0]->qual) )
        {
            if ( bcf_float_is_missing(out->qual) || out->qual < files->readers[i].buffer[0]->qual ) out->qual = files->readers[i].buffer[0]->qual;
        }
    }
    //QUAL is ignored in GATK GenotypeGVCF
    if(g_do_gatk_merge)
        bcf_float_set_missing(out->qual);

    // set ID
    if ( !tmps->l ) kputs(".", tmps);
    if ( out->d.id ) free(out->d.id);
    out->d.id = strdup(tmps->s);

    // set alleles
    ma->nout_als = 0;
    for (i=1; i<ma->nals; i++)  //0 is reference
    {
        if ( !al_idxs[i] ) continue;
        ma->nout_als++;

        // Adjust the indexes, the allele map could be created for multiple collapsed records,
        //  some of which might be unused for this output line
        int ir, j;
        for (ir=0; ir<files->nreaders; ir++)
        {
            if ( !ma->has_line[ir] ) continue;
            bcf1_t *line = files->readers[ir].buffer[0];
            for (j=1; j<line->n_allele; j++)    //0 is reference
                if ( ma->d[ir][0].map[j]==i ) ma->d[ir][0].map[j] = ma->nout_als;
        }
    }
    // Expand the arrays and realloc the alleles string. Note that all alleles are in a single allocated block.
    ma->nout_als++;
    hts_expand0(char*, ma->nout_als, ma->mout_als, ma->out_als);
    int k = 0;
    for (i=0; i<ma->nals; i++)
        if ( i==0 || al_idxs[i] ) ma->out_als[k++] = strdup(ma->als[i]);
    assert( k==ma->nout_als );
    normalize_alleles(ma->out_als, ma->nout_als);
    bcf_update_alleles(out_hdr, out, (const char**) ma->out_als, ma->nout_als);
#ifdef USE_ID_MAP
    if(out->n_allele == 2 && out->d.var[1].type & VCF_NON_REF)
        args->m_idmap.m_merged_has_only_non_ref = 1;
    else
        args->m_idmap.m_merged_has_only_non_ref = 0;
#endif
    free(al_idxs);
    for (i=0; i<ma->nout_als; i++) free(ma->out_als[i]);
}

void merge_filter(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->d.n_flt = 0;
    if(!g_do_gatk_merge) //GATK GenotypeGVCF merge ignores filters
        for (i=0; i<files->nreaders; i++)
        {
            if ( !ma->has_line[i]) continue;

            bcf_sr_t *reader = &files->readers[i];
            bcf1_t *line = reader->buffer[0];
            bcf_hdr_t *hdr = reader->processed_header;
            bcf_unpack(line, BCF_UN_ALL);

            int k;
            for (k=0; k<line->d.n_flt; k++)
            {
                const char *flt = hdr->id[BCF_DT_ID][line->d.flt[k]].key;
                kitr = kh_get(strdict, tmph, flt);
                //if filter not already added, add
                if ( kitr == kh_end(tmph) )
                {
                    int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, flt);
                    if ( id==-1 ) error("The filter not defined: %s\n", flt);
                    hts_expand(int,out->d.n_flt+1,ma->mflt,ma->flt);
                    ma->flt[out->d.n_flt] = id;
                    out->d.n_flt++;
                    kh_put(strdict, tmph, flt, &ret);
                }
            }
        }
    // Check if PASS is not mixed with other filters
    if ( out->d.n_flt>1 )
    {
        int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "PASS");
        for (i=0; i<out->d.n_flt; i++)
            if ( ma->flt[i]==id ) break;
        if ( i<out->d.n_flt )
        {
            out->d.n_flt--;
            for (; i<out->d.n_flt; i++) ma->flt[i] = ma->flt[i+1];
        }
    }
    out->d.flt = ma->flt;
}

static void bcf_info_set_id(bcf1_t *line, bcf_info_t *info, int id, kstring_t *tmp_str)
{
    assert( !info->vptr_free );

    uint8_t *ptr = info->vptr - info->vptr_off;
    bcf_dec_typed_int1(ptr, &ptr);

    tmp_str->l = 0;
    bcf_enc_int1(tmp_str, id);

    if ( tmp_str->l == ptr - info->vptr + info->vptr_off )
    {
        // the new id is represented with the same number of bytes
        memcpy(info->vptr - info->vptr_off, tmp_str->s, tmp_str->l);
        return;
    }

    kputsn_(ptr, info->vptr - ptr, tmp_str);
    info->vptr_off = tmp_str->l;
    kputsn_(info->vptr, info->len << bcf_type_shift[info->type], tmp_str);

    info->vptr = (uint8_t*) tmp_str->s + info->vptr_off;
    info->vptr_free = 1;
    line->d.shared_dirty |= BCF1_DIRTY_INF;
    tmp_str->s = NULL;
    tmp_str->m = 0;
    tmp_str->l = 0;
}

int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0;    // i-th field in src string
    while ( ith_src<isrc && start_src<src_len )
    {
        if ( src[start_src]==',' ) { ith_src++; }
        start_src++;
    }
    if ( ith_src!=isrc ) return -1; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' ) return 0;   // don't write missing values, dst is already initialized

    int ith_dst = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l )
    {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    if ( ith_dst!=idst ) return -2;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_dst]!=',' ) end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' ) return 0;   // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move  = dst->l - end_dst + 1;  // how many bytes must be moved (including \0)
    if ( ndst_shift )
    {
        ks_resize(dst, dst->l + ndst_shift + 1);    // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
    return 0;
}

static void merge_AGR_info_tag(bcf_hdr_t *hdr, bcf1_t *line, bcf_info_t *info, int len, maux1_t *als, AGR_info_t *agr)
{
    int i;
    maux1_t* als = &(ma->d[reader_idx][0]);
    if ( !agr->nbuf )
    {
        if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
        {
            agr->nbuf = agr->nvals*sizeof(int); //4==sizeof(int)==sizeof(float)
            hts_expand(uint8_t,agr->nbuf,agr->mbuf,agr->buf);
            if ( info->type!=BCF_BT_FLOAT )
            {
                int32_t *tmp = (int32_t*) agr->buf;
                for (i=0; i<agr->nvals; i++) tmp[i] = bcf_int32_missing;
            }
            else
            {
                float *tmp = (float*) agr->buf;
                for (i=0; i<agr->nvals; i++) bcf_float_set_missing(tmp[i]);
            }
        }
        else if ( info->type==BCF_BT_CHAR )
        {
            kstring_t tmp; tmp.l = 0; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
            kputc('.',&tmp);
            for (i=1; i<agr->nvals; i++) kputs(",.",&tmp);
            agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
        }
        else
            error("Not ready for type [%d]: %s at %d\n", info->type,agr->hdr_tag,line->pos+1);
    }

    if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
    {
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int ifrom = len==BCF_VL_A ? 1 : 0;  //0 is REF allele for BCF_VL_A
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori, inew; \
                if(line->m_non_ref_idx >= 0) \
                { \
                    initialize_ARG_vectors_for_NON_REF(line, len, als, (ma->nout_als), src, (agr->nvals), tgt, out_type_t); \
                } \
                for (iori=ifrom; iori<line->n_allele; iori++) \
                { \
                    if ( is_vector_end ) break; \
                    if ( is_missing ) continue; \
                    inew = als->map[iori] - ifrom; \
                    tgt[inew] = *src; \
                    src++; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  *src==bcf_int8_missing, *src==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
        else    //G type
        {
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori,jori, inew,jnew; \
                if(line->m_non_ref_idx >= 0)  /*NON_REF allele*/ \
                { \
                    initialize_ARG_vectors_for_NON_REF(line, len, als, (ma->nout_als), src, (agr->nvals), tgt, out_type_t); \
                } \
                for (iori=0; iori<line->n_allele; iori++) \
                { \
                    inew = als->map[iori]; \
                    for (jori=0; jori<=iori; jori++) \
                    { \
                        jnew = als->map[jori]; \
                        int kori = iori*(iori+1)/2 + jori; \
                        if ( is_vector_end ) break; \
                        if ( is_missing ) continue; \
                        int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                        tgt[knew] = src[kori]; \
                    } \
                    if ( jori<=iori ) break; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  src[kori]==bcf_int8_missing,  src[kori]==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, src[kori]==bcf_int16_missing, src[kori]==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, src[kori]==bcf_int32_missing, src[kori]==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(src[kori]), bcf_float_is_vector_end(src[kori]), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
    }
    else
    {
        kstring_t tmp; tmp.l = agr->nbuf; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int iori, ifrom = len==BCF_VL_A ? 1 : 0;
            for (iori=ifrom; iori<line->n_allele; iori++)
            {
                int ret = copy_string_field((char*)info->vptr, iori-ifrom, info->len, &tmp, als->map[iori]-ifrom);
                if ( ret )
                    error("Error at %s:%d: wrong number of fields in %s?\n", bcf_seqname(hdr,line),line->pos+1,agr->hdr_tag);
            }
        }
        else
        {
            int iori,jori, inew,jnew;
            for (iori=0; iori<line->n_allele; iori++)
            {
                inew = als->map[iori];
                for (jori=0; jori<=iori; jori++)
                {
                    jnew = als->map[jori];
                    int kori = iori*(iori+1)/2 + jori;
                    int knew = bcf_alleles2gt(inew,jnew);
                    int ret  = copy_string_field((char*)info->vptr, kori, info->len, &tmp, knew);
                    if ( ret )
                        error("Error at %s:%d: wrong number of fields in %s?\n", bcf_seqname(hdr,line),line->pos+1,agr->hdr_tag);
                }
            }
        }
        agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
    }
}

void merge_info(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, j;
#ifndef USE_ID_MAP
    int ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
#else
    args->m_idmap.m_merged_DP_info_idx = -1;
    memset(args->m_idmap.m_DP_info_vals, -1, args->files->nreaders*sizeof(int));
#endif
    maux_t *ma = args->maux;
    ma->nAGR_info = 0;
    out->n_info   = 0;
    info_rules_reset(args);     //zero out everything
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->processed_header;
#ifdef USE_ID_MAP
        reader_idmap* curr_map = &(args->m_idmap.m_readers_map[i]);
#endif
        for (j=0; j<line->n_info; j++) 
        {
            bcf_info_t *inf = &line->d.info[j];
            if(inf->vptr == 0)
                continue;

            const char *key = hdr->id[BCF_DT_ID][inf->key].key;
            if ( !strcmp("AC",key) || !strcmp("AN",key) ) continue;  // AC and AN are done in merge_format() after genotypes are done
#ifdef USE_ID_MAP
            ASSERT(inf->key < curr_map->m_num_ids);
            int id = curr_map->m_id_2_merged_id[inf->key];
            ASSERT(id >= 0 && id < out_hdr->n[BCF_DT_ID] && id < args->m_idmap.m_num_merged_ids);
            int* merged_idx_ptr = &(args->m_idmap.m_merged_id_2_merged_idx[id]);
#ifdef DEBUG
            ASSERT(bcf_hdr_idinfo_exists(out_hdr, BCF_HL_INFO, id));    //should be INFO field 
            const char* out_key = bcf_hdr_int2id(out_hdr, BCF_DT_ID, id);
            ASSERT(out_key && strcmp(out_key, key) == 0);       //same as this key
#endif //DEBUG
            //is END info tag, but pos == end - don't add END tag
            if(id == args->m_idmap.m_merged_END_info_id && bcf_get_end_point(line) == out->pos)
                continue;
#else  //USE_ID_MAP
            int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, key);
            if ( id==-1 ) error("Error: The INFO field is not defined in the header: %s\n", key);

            kitr = kh_get(strdict, tmph, key);  // have we seen the tag in one of the readers?
#endif //USE_ID_MAP 
            int len = bcf_hdr_id2length(hdr,BCF_HL_INFO,inf->key);
            if ( args->nrules )
            {
#ifdef USE_ID_MAP
                //avoid string search - use indexing
                int info_rule_idx = args->m_idmap.m_merged_id_2_info_rule_idx[id];
                ASSERT(info_rule_idx < args->nrules);
                info_rule_t* rule = (info_rule_idx < 0) ? 0 : &(args->rules[info_rule_idx]);
#else
                info_rule_t *rule = (info_rule_t*) bsearch(key, args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key);
#endif
                if ( rule ) 
                {
                    ASSERT(strcmp(rule->hdr_tag, out_key) == 0);
                    maux1_t *als = ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R ) ? &ma->d[i][0] : NULL;
                    if ( info_rules_add_values(args, hdr, line, rule, als, len) )
                    {
#ifdef USE_ID_MAP
                        if(id == args->m_idmap.m_merged_DP_info_id)     //is DP info field
                            args->m_idmap.m_DP_info_vals[i] = ((int32_t*)(rule->vals))[rule->nvals-1];
#endif
                        continue;
                    }
                }
            }

            // TODO: Number=AGR tags should use the newer info_rules_* functions (info_rules_merge_first to be added)
            // and merge_AGR_info_tag to be made obsolete.
            if ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R  ) // Number=R,G,A requires special treatment
            {
#ifdef USE_ID_MAP
                if(*merged_idx_ptr == -1)
#else
                if ( kitr == kh_end(tmph) )
#endif
                {
                    // first occurance in this reader, alloc arrays
                    ma->nAGR_info++;
                    hts_expand0(AGR_info_t,ma->nAGR_info,ma->mAGR_info,ma->AGR_info);
#ifdef USE_ID_MAP
                    *merged_idx_ptr = ma->nAGR_info - 1;
#else
                    kitr = kh_put(strdict, tmph, key, &ret);
                    kh_val(tmph,kitr) = ma->nAGR_info - 1;
#endif
                    ma->AGR_info[ma->nAGR_info-1].hdr_tag = key;
                    ma->AGR_info[ma->nAGR_info-1].type  = bcf_hdr_id2type(hdr,BCF_HL_INFO,inf->key);
                    ma->AGR_info[ma->nAGR_info-1].nbuf  = 0;    // size of the buffer
                    switch (len)
                    {
                        case BCF_VL_A: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als - 1; break;
                        case BCF_VL_G: ma->AGR_info[ma->nAGR_info-1].nvals = bcf_alleles2gt(ma->nout_als-1,ma->nout_als-1)+1; break;
                        case BCF_VL_R: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als; break;
                    }
                }
#ifdef USE_ID_MAP
                int idx = *merged_idx_ptr;
#ifdef DEBUG
                ASSERT(idx >= 0 && idx < ma->nAGR_info);
                ASSERT(strcmp(ma->AGR_info[idx].hdr_tag, key) == 0);    //same as key
#endif
#else   //USE_ID_MAP
                kitr = kh_get(strdict, tmph, key);
                int idx = kh_val(tmph, kitr);
                if ( idx<0 ) error("Error occurred while processing INFO tag \"%s\" at %s:%d\n", key,bcf_seqname(hdr,line),line->pos+1);
#endif  //USE_ID_MAP
                merge_AGR_info_tag(hdr, line,inf,len,&ma->d[i][0],&ma->AGR_info[idx]);
                continue;
            }
#ifdef USE_ID_MAP
            if(*merged_idx_ptr == -1)
#else
            if ( kitr == kh_end(tmph) )
#endif
            {
                hts_expand0(bcf_info_t,out->n_info+1,ma->minf,ma->inf);
                ma->inf[out->n_info].key  = id;
                ma->inf[out->n_info].type = inf->type;
                ma->inf[out->n_info].len  = inf->len;
                ma->inf[out->n_info].vptr = inf->vptr;
                ma->inf[out->n_info].v1.i = inf->v1.i;
                ma->inf[out->n_info].v1.f = inf->v1.f;
                ma->inf[out->n_info].vptr_off  = inf->vptr_off;
                ma->inf[out->n_info].vptr_len  = inf->vptr_len;
                ma->inf[out->n_info].vptr_free = inf->vptr_free;
                if ( (args->output_type & FT_BCF) && id!=bcf_hdr_id2int(hdr, BCF_DT_ID, key) )
                {
                    // The existing packed info cannot be reused. Change the id.
                    // Although quite hacky, it's faster than anything else given
                    // the data structures
                    bcf_info_set_id(out, &ma->inf[out->n_info], id, &args->tmps);
                }
                out->n_info++;
#ifdef USE_ID_MAP
                *merged_idx_ptr = (out->n_info - 1);
#else
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_val(tmph,kitr) = -(out->n_info-1);   // arbitrary negative value TODO:WHY??
#endif
            }
        }
    }
    out->d.info = ma->inf;
    out->d.m_info = ma->minf;
    for (i=0; i<args->nrules; i++)
    {
        args->rules[i].merger(args->out_hdr, out, &args->rules[i]);
        //DP field was added by merger
        if(args->rules[i].nvals && args->rules[i].merged_id == args->m_idmap.m_merged_DP_info_id)
        {
            args->m_idmap.m_merged_DP_info_idx = out->n_info-1;//last added
            ASSERT(strcmp(bcf_hdr_int2id(out_hdr, BCF_DT_ID, out->d.info[args->m_idmap.m_merged_DP_info_idx].key),"DP") == 0);
        }
    }
    for (i=0; i<ma->nAGR_info; i++)
    {
        AGR_info_t *agr = &ma->AGR_info[i];
        bcf_update_info(out_hdr,out,agr->hdr_tag,agr->buf,agr->nvals,agr->type);
    }
}

void update_AN_AC(bcf_hdr_t *hdr, bcf1_t *line)
{
    int32_t an = 0, *tmp = (int32_t*) malloc(sizeof(int)*line->n_allele);
    int ret = bcf_calc_ac(hdr, line, tmp, BCF_UN_FMT);
    if ( ret>0 )
    {
        int i;
        for (i=0; i<line->n_allele; i++) an += tmp[i];
        bcf_update_info_int32(hdr, line, "AN", &an, 1);
        bcf_update_info_int32(hdr, line, "AC", tmp+1, line->n_allele-1);
    }
    free(tmp);
}

void merge_GT(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);
   
    //nsize - max entries per sample (max ploidy here)
    //see definition of fmt_map for why nsize is ploidy
    int max_ploidy = 0, msize = sizeof(int32_t);
    for (i=0; i<files->nreaders; i++)
    {
        if ( !fmt_map[i] ) continue;
        if ( fmt_map[i]->n > max_ploidy ) max_ploidy = fmt_map[i]->n;
    }

    if ( ma->ntmp_arr < nsamples*max_ploidy*msize )
    {
        ma->ntmp_arr = nsamples*max_ploidy*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }
    memset(ma->smpl_ploidy,0,nsamples*sizeof(int));

    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->processed_header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        int32_t *tmp  = (int32_t *) ma->tmp_arr + ismpl*max_ploidy;

        int j, k;
        if ( !fmt_ori )
        {
            // missing values: assume maximum ploidy
            for (j=0; j<bcf_hdr_nsamples(hdr); j++)
            {
                for (k=0; k<max_ploidy; k++) { tmp[k] = 0; ma->smpl_ploidy[ismpl+j]++; }
                tmp += max_ploidy;
            }
            ismpl += bcf_hdr_nsamples(hdr);
            continue;
        }

        #define BRANCH(type_t, vector_end) { \
            type_t *p_ori  = (type_t*) fmt_ori->p; \
            if ( !ma->d[i][0].als_differ ) \
            { \
                /* the allele numbering is unchanged */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (k=0; k<fmt_ori->n; k++) \
                    { \
                        if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                        ma->smpl_ploidy[ismpl+j]++; \
                        if ( bcf_gt_is_missing(p_ori[k]) ) tmp[k] = 0; /* missing allele */ \
                        else tmp[k] = p_ori[k]; \
                    } \
                    for (; k<max_ploidy; k++) tmp[k] = bcf_int32_vector_end; \
                    tmp += max_ploidy; \
                    p_ori += fmt_ori->n; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
            { \
                for (k=0; k<fmt_ori->n; k++) \
                { \
                    if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                    ma->smpl_ploidy[ismpl+j]++; \
                    if ( bcf_gt_is_missing(p_ori[k]) ) tmp[k] = 0; /* missing allele */ \
                    else \
                    { /*WTF??*/ \
                        int al = (p_ori[k]>>1) - 1; \
                        al = al<=0 ? al + 1 : ma->d[i][0].map[al] + 1; \
                        tmp[k] = (al << 1) | ((p_ori[k])&1); \
                    } \
                } \
                for (; k<max_ploidy; k++) tmp[k] = bcf_int32_vector_end; \
                tmp += max_ploidy; \
                p_ori += fmt_ori->n; \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (fmt_ori->type)
        {
            case BCF_BT_INT8: BRANCH(int8_t,   bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
            default: error("Unexpected case: %d\n", fmt_ori->type);
        }
        #undef BRANCH
    }
    bcf_update_format_int32(out_hdr, out, "GT", (int32_t*)ma->tmp_arr, nsamples*max_ploidy);
}

void merge_format_field(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);

    const char *key = NULL;
    int nsize = 0, length = BCF_VL_FIXED, type = -1;

    int out_id = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        if ( !fmt_map[i] ) continue;
        if ( !key ) key = files->readers[i].processed_header->id[BCF_DT_ID][fmt_map[i]->id].key;
        type = fmt_map[i]->type;
        int length_descriptor = bcf_hdr_id2length(files->readers[i].processed_header, BCF_HL_FMT, fmt_map[i]->id);
#ifdef USE_ID_MAP
        out_id = args->m_idmap.m_readers_map[i].m_id_2_merged_id[fmt_map[i]->id];
#endif
        if ( IS_VL_G(files->readers[i].processed_header, fmt_map[i]->id) )
        { 
            length = BCF_VL_G; 
            nsize = out->n_allele*(out->n_allele + 1)/2; 
            break;
        }
        if (length_descriptor == BCF_VL_A || length_descriptor == BCF_VL_R)
        {
            length = length_descriptor;
            nsize = (length_descriptor == BCF_VL_A) ? out->n_allele - 1 : out->n_allele;
            break;
        }
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    int msize = sizeof(float)>sizeof(int32_t) ? sizeof(float) : sizeof(int32_t);
    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }

    // Fill the temp array for all samples by collecting values from all files
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->processed_header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        maux1_t* als = &(ma->d[i][0]);
        if ( fmt_ori )
        {
            type = fmt_ori->type;
            int nals_ori = reader->buffer[0]->n_allele;
            if ( length==BCF_VL_G )
            {
                // if all fields are missing then n==1 is valid
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori*(nals_ori+1)/2 && fmt_map[i]->n != nals_ori )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
            else if ( length==BCF_VL_A )
            {
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori-1 )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
            else if ( length==BCF_VL_R )
            {
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
        }

        // set the values
        #define BRANCH(tgt_type_t, src_type_t, src_is_missing, src_is_vector_end, tgt_set_missing, tgt_set_vector_end) { \
            int j, l, k; \
            tgt_type_t *tgt = (tgt_type_t *) ma->tmp_arr + ismpl*nsize; \
            if ( !fmt_ori ) \
            { \
                /* the field is not present in this file, set missing values */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt_set_missing; tgt++; for (l=1; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            assert( ma->has_line[i] ); \
            bcf1_t *line    = reader->buffer[0]; \
            src_type_t *src = (src_type_t*) fmt_ori->p; \
            if ( (length!=BCF_VL_G && length!=BCF_VL_A && length!=BCF_VL_R) || (line->n_allele==out->n_allele && !ma->d[i][0].als_differ) ) \
            { \
                /* alleles unchanged, copy over */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (l=0; l<fmt_ori->n; l++) \
                    { \
                        if ( src_is_vector_end ) break; \
                        else if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        tgt++; src++; \
                    } \
                    for (k=l; k<nsize; k++) { tgt_set_vector_end; tgt++; } \
                    src += fmt_ori->n - l; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            if ( length==BCF_VL_G ) \
            { \
                /* Number=G tags */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    src = (src_type_t*) fmt_ori->p + j*fmt_ori->n; \
                    if ( (src_is_missing && fmt_ori->n==1) || (++src && src_is_vector_end) ) \
                    { \
                        /* tag with missing value "." */ \
                        tgt_set_missing; \
                        for (l=1; l<nsize; l++) { tgt++; tgt_set_vector_end; } \
                        continue; \
                    } \
                    int ngsize = ma->smpl_ploidy[ismpl+j]==1 ? out->n_allele : out->n_allele*(out->n_allele + 1)/2; \
                    /*printf("Out num alleles %d num GT %d input num allele %d non_ref idx %d\n",out->n_allele, ngsize, line->n_allele, line->m_non_ref_idx);*/ \
                    if(line->m_non_ref_idx >= 0)  /*NON_REF allele*/ \
                    { \
                        int src_length = (ma->smpl_ploidy[ismpl+j]==1) ? line->n_allele : (line->n_allele*(line->n_allele+1))/2; \
                        initialize_ARG_vectors_for_NON_REF(line, length, als, (out->n_allele), src, src_length, tgt, tgt_type_t); \
                    } \
                    else \
                    { \
                        for (l=0; l<ngsize; l++) { tgt_set_missing; tgt++; } \
                        for (; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                    } \
                    if ( ma->smpl_ploidy[ismpl+j]==1 ) \
                    { \
                        /* Haploid */ \
                        int iori, inew; \
                        for (iori=0; iori<line->n_allele; iori++) \
                        { \
                            inew = ma->d[i][0].map[iori]; \
                            src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + iori; \
                            tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                            if ( src_is_vector_end ) break; \
                            if ( src_is_missing ) tgt_set_missing; \
                            else *tgt = *src; \
                        } \
                    } \
                    else \
                    { \
                        /* Diploid */ \
                        int iori,jori, inew,jnew; \
                        for (iori=0; iori<line->n_allele; iori++) \
                        { \
                            inew = ma->d[i][0].map[iori]; \
                            for (jori=0; jori<=iori; jori++) \
                            { \
                                jnew = ma->d[i][0].map[jori]; \
                                int kori = iori*(iori+1)/2 + jori; \
                                int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                                src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + kori; \
                                tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + knew; \
                                if ( src_is_vector_end ) \
                                { \
                                    iori = line->n_allele; \
                                    break; \
                                } \
                                if ( src_is_missing ) tgt_set_missing; \
                                else *tgt = *src; \
                            } \
                        } \
                    } \
                } \
            } \
            else \
            { \
                /* Number=A or Number=R tags */ \
                int ifrom = (length == BCF_VL_A) ? 1 : 0; \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    src = (src_type_t*) (fmt_ori->p + j*fmt_ori->size); \
                    if ( (src_is_missing && fmt_ori->n==1) || (++src && src_is_vector_end) ) \
                    { \
                        /* tag with missing value "." */ \
                        tgt_set_missing; \
                        for (l=1; l<nsize; l++) { tgt++; tgt_set_vector_end; } \
                        continue; \
                    } \
                    src = (src_type_t*) (fmt_ori->p + j*fmt_ori->size); \
                    if(line->m_non_ref_idx >= 0)  /*NON_REF allele*/ \
                    { \
                        int src_length =  (length == BCF_VL_R) ? line->n_allele : line->n_allele-1; \
                        initialize_ARG_vectors_for_NON_REF(line, length, als, (out->n_allele), src, src_length, tgt, tgt_type_t); \
                    } \
                    else \
                        for (l=0; l<nsize; l++) { tgt_set_missing; tgt++; } \
                    int iori,inew; \
                    for (iori=ifrom; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori] - ifrom; \
                        tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                        if ( src_is_vector_end ) break; \
                        if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        src++; \
                    } \
                } \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (type)
        {
            case BCF_BT_INT8:  BRANCH(int32_t,  int8_t, *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT16: BRANCH(int32_t, int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_FLOAT: BRANCH(float, float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*tgt), bcf_float_set_vector_end(*tgt)); break;
            case BCF_BT_CHAR:  BRANCH(uint8_t, uint8_t, *src==bcf_str_missing, *src==bcf_str_vector_end, *tgt=bcf_str_missing, *tgt=bcf_str_vector_end); break;
            default: error("Unexpected case: %d, %s\n", type, key);
        }
        #undef BRANCH
    }
    if ( type==BCF_BT_FLOAT )
        bcf_update_format_float(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else if ( type==BCF_BT_CHAR )
        bcf_update_format_char(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else
        bcf_update_format_int32(out_hdr, out, key, (int32_t*)ma->tmp_arr, nsamples*nsize);
}

#ifdef USE_ID_MAP
void update_DP_info(args_t* args, bcf1_t* out)
{
    //DP was added as an INFO field in this line
    if(g_do_gatk_merge && args->m_idmap.m_merged_DP_info_idx >= 0)
    {
        int32_t* info_ptr = args->m_idmap.m_DP_info_vals;
        int32_t sum = 0;
        int i = 0;
        int j = 0;
        int k = 0;

#define BRANCH(MIN_DP_type_t, DP_type_t, is_missing_MIN_DP, is_terminator_MIN_DP, is_missing_DP, is_terminator_DP) \
        { \
            MIN_DP_type_t* format_MIN_DP_ptr = (args->m_idmap.m_merged_MIN_DP_format_idx >= 0) ? \
                (MIN_DP_type_t*)(out->d.fmt[args->m_idmap.m_merged_MIN_DP_format_idx].p) : 0; \
            DP_type_t* format_DP_ptr = (args->m_idmap.m_merged_DP_format_idx >= 0) ? \
                (DP_type_t*)(out->d.fmt[args->m_idmap.m_merged_DP_format_idx].p) : 0; \
            for(i=0;i<args->files->nreaders;++i) \
            { \
                int local_nsamples = bcf_hdr_nsamples(args->files->readers[i].processed_header); \
                if(info_ptr[i] < 0)    /*INFO field DP missing, use value from FORMAT field*/ \
                { \
                    int32_t local_sum = 0; /*Compute mean over DP values in FORMAT fields*/ \
                    for(j=0;j<local_nsamples;++j) \
                    { \
                        if(format_MIN_DP_ptr == 0 || is_missing_MIN_DP || is_terminator_MIN_DP) \
                            if(format_DP_ptr == 0 || is_missing_DP || is_terminator_DP) \
                                continue; \
                            else \
                                local_sum += (int32_t)(format_DP_ptr[k+j]); \
                        else \
                            local_sum += (int32_t)(format_MIN_DP_ptr[k+j]); \
                    } \
                    sum += (local_sum/local_nsamples); \
                } \
                else \
                    sum += info_ptr[i]; \
                k += local_nsamples; \
            } \
        }
        int MIN_DP_type = (args->m_idmap.m_merged_MIN_DP_format_idx >= 0) ?
            out->d.fmt[args->m_idmap.m_merged_MIN_DP_format_idx].type : BCF_BT_INT32;
        int DP_type = (args->m_idmap.m_merged_DP_format_idx >= 0) ?
            out->d.fmt[args->m_idmap.m_merged_DP_format_idx].type : BCF_BT_INT32;
        //need to consider all combinations of MIN_DP_type, DP_type - this is as good as any
        switch(MIN_DP_type << 16 | DP_type)
        {
            case ((BCF_BT_INT8 << 16) | BCF_BT_INT8):
                BRANCH(int8_t, int8_t, (format_MIN_DP_ptr[k+j] == bcf_int8_missing), (format_MIN_DP_ptr[k+j] == bcf_int8_vector_end), (format_DP_ptr[k+j] == bcf_int8_missing), (format_DP_ptr[k+j] == bcf_int8_vector_end));
                break;
            case ((BCF_BT_INT8 << 16) | BCF_BT_INT16):
                BRANCH(int8_t, int16_t, (format_MIN_DP_ptr[k+j] == bcf_int8_missing), (format_MIN_DP_ptr[k+j] == bcf_int8_vector_end), (format_DP_ptr[k+j] == bcf_int16_missing), (format_DP_ptr[k+j] == bcf_int16_vector_end));
                break;
            case ((BCF_BT_INT8 << 16) | BCF_BT_INT32):
                BRANCH(int8_t, int32_t, (format_MIN_DP_ptr[k+j] == bcf_int8_missing), (format_MIN_DP_ptr[k+j] == bcf_int8_vector_end), (format_DP_ptr[k+j] == bcf_int32_missing), (format_DP_ptr[k+j] == bcf_int32_vector_end));
                break;
            case ((BCF_BT_INT16 << 16) | BCF_BT_INT8):
                BRANCH(int16_t, int8_t, (format_MIN_DP_ptr[k+j] == bcf_int16_missing), (format_MIN_DP_ptr[k+j] == bcf_int16_vector_end), (format_DP_ptr[k+j] == bcf_int8_missing), (format_DP_ptr[k+j] == bcf_int8_vector_end));
                break;
            case ((BCF_BT_INT16 << 16) | BCF_BT_INT16):
                BRANCH(int16_t, int16_t, (format_MIN_DP_ptr[k+j] == bcf_int16_missing), (format_MIN_DP_ptr[k+j] == bcf_int16_vector_end), (format_DP_ptr[k+j] == bcf_int16_missing), (format_DP_ptr[k+j] == bcf_int16_vector_end));
                break;
            case ((BCF_BT_INT16 << 16) | BCF_BT_INT32):
                BRANCH(int16_t, int32_t, (format_MIN_DP_ptr[k+j] == bcf_int16_missing), (format_MIN_DP_ptr[k+j] == bcf_int16_vector_end), (format_DP_ptr[k+j] == bcf_int32_missing), (format_DP_ptr[k+j] == bcf_int32_vector_end));
                break;
            case ((BCF_BT_INT32 << 16) | BCF_BT_INT8):
                BRANCH(int32_t, int8_t, (format_MIN_DP_ptr[k+j] == bcf_int32_missing), (format_MIN_DP_ptr[k+j] == bcf_int32_vector_end), (format_DP_ptr[k+j] == bcf_int8_missing), (format_DP_ptr[k+j] == bcf_int8_vector_end));
                break;
            case ((BCF_BT_INT32 << 16) | BCF_BT_INT16):
                BRANCH(int32_t, int16_t, (format_MIN_DP_ptr[k+j] == bcf_int32_missing), (format_MIN_DP_ptr[k+j] == bcf_int32_vector_end), (format_DP_ptr[k+j] == bcf_int16_missing), (format_DP_ptr[k+j] == bcf_int16_vector_end));
                break;
            case ((BCF_BT_INT32 << 16) | BCF_BT_INT32):
                BRANCH(int32_t, int32_t, (format_MIN_DP_ptr[k+j] == bcf_int32_missing), (format_MIN_DP_ptr[k+j] == bcf_int32_vector_end), (format_DP_ptr[k+j] == bcf_int32_missing), (format_DP_ptr[k+j] == bcf_int32_vector_end));
                break;
            default:
                assert(0);
        }
#undef BRANCH
        ASSERT(args->m_idmap.m_merged_DP_info_idx >= 0 && args->m_idmap.m_merged_DP_info_idx < out->n_info);
        out->d.info[args->m_idmap.m_merged_DP_info_idx].v1.i = sum;
    }
}
#endif

void merge_format(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    //map from ith format entry of merged line to jth sample's corresponding bcf_fmt_t* entry 
    //ma->fmt_map[i*nreaders+j] gives the bcf_fmt_t* entry for sample j 
    if ( !ma->nfmt_map ) 
    {
        ma->nfmt_map = 2;
        ma->fmt_map  = (bcf_fmt_t**) calloc(ma->nfmt_map*files->nreaders, sizeof(bcf_fmt_t*));
    }
    else
        memset(ma->fmt_map, 0, ma->nfmt_map*files->nreaders*sizeof(bcf_fmt_t**));

#ifndef USE_ID_MAP
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    int ret = 0;
#else
    args->m_idmap.m_merged_DP_format_idx = -1;
    args->m_idmap.m_merged_MIN_DP_format_idx = -1;
#endif
    int i, j, has_GT = 0, max_ifmt = 0; // max fmt index

    //Get all format fields and collect info about allele re-numbering for every sample
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->processed_header;
#ifdef USE_ID_MAP
        reader_idmap* curr_map = &(args->m_idmap.m_readers_map[i]);
#endif
        for (j=0; j<line->n_fmt; j++) 
        {
            // Wat this tag already seen?
            bcf_fmt_t *fmt = &line->d.fmt[j];
            if(fmt->p == 0)     //deleted
                continue;
            const char *key = hdr->id[BCF_DT_ID][fmt->id].key;
#ifdef USE_ID_MAP
            //For lines with variants other than NON_REF, GQ format field is dropped
            if(g_do_gatk_merge && !(args->m_idmap.m_merged_has_only_non_ref) && fmt->id == curr_map->m_GQ_id)
                continue;
            ASSERT(fmt->id < curr_map->m_num_ids);
            int out_id = curr_map->m_id_2_merged_id[fmt->id];
            ASSERT(out_id >= 0 && out_id < out_hdr->n[BCF_DT_ID] && out_id < args->m_idmap.m_num_merged_ids);
            int* merged_idx_ptr = &(args->m_idmap.m_merged_id_2_merged_idx[out_id]);
            int ifmt = *merged_idx_ptr;
#ifdef DEBUG
            ASSERT(bcf_hdr_idinfo_exists(out_hdr, BCF_HL_FMT, out_id));    //should be FORMAT field 
            const char* out_key = bcf_hdr_int2id(out_hdr, BCF_DT_ID, out_id);
            ASSERT(out_key && strcmp(out_key, key) == 0);       //same as this key
#endif
#else   //USE_ID_MAP
            kitr = kh_get(strdict, tmph, key);
            int ifmt = -1;
            if ( kitr != kh_end(tmph) )
                ifmt = kh_value(tmph, kitr);    // seen
#endif  //USE_ID_MAP
            if(ifmt == -1)
            {
                // new FORMAT tag
                if ( key[0]=='G' && key[1]=='T' && key[2]==0 ) { has_GT = 1; ifmt = 0; }
                else
                {
                    ifmt = ++max_ifmt;  //GT is fmt index 0
#ifdef USE_ID_MAP
                    if(out_id == args->m_idmap.m_merged_DP_format_id)
                        args->m_idmap.m_merged_DP_format_idx = ifmt;
                    else
                        if(out_id == args->m_idmap.m_merged_MIN_DP_format_id)
                            args->m_idmap.m_merged_MIN_DP_format_idx = ifmt;
#endif
                    if ( max_ifmt >= ma->nfmt_map )
                    {
                        ma->fmt_map = (bcf_fmt_t**) realloc(ma->fmt_map, sizeof(bcf_fmt_t*)*(max_ifmt+1)*files->nreaders);
                        memset(ma->fmt_map+ma->nfmt_map*files->nreaders, 0, (max_ifmt-ma->nfmt_map+1)*files->nreaders*sizeof(bcf_fmt_t*));
                        ma->nfmt_map = max_ifmt+1;
                    }
                }
#ifdef USE_ID_MAP
                *merged_idx_ptr = ifmt;
#else
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_value(tmph, kitr) = ifmt;
#endif
            }
            ma->fmt_map[ifmt*files->nreaders+i] = fmt;
        }
        // Check if the allele numbering must be changed
        for (j=1; j<reader->buffer[0]->n_allele; j++)   //0 is reference
            if ( ma->d[i][0].map[j]!=j ) break;
        ma->d[i][0].als_differ = j==reader->buffer[0]->n_allele ? 0 : 1;
    }

    out->n_sample = bcf_hdr_nsamples(out_hdr);
    if ( has_GT )
        merge_GT(args, ma->fmt_map, out);
    update_AN_AC(out_hdr, out);

    //0 is GT field - already merged
    for (i=1; i<=max_ifmt; i++)
        merge_format_field(args, &ma->fmt_map[i*files->nreaders], out);
    //bcftools and htslib clusterf***
    if ( out->d.info!=ma->inf )
    {
        // hacky, we rely on htslib internals: bcf_update_info() reallocated the info
        ma->inf  = out->d.info;
        ma->minf = out->d.m_info;
    }
    update_DP_info(args, out);

    out->d.indiv_dirty = 1;
}

// The core merging function, one or none line from each reader
void merge_line(args_t *args)
{
    bcf1_t *out = args->out_line;
    bcf_clear1(out);
    out->unpacked = BCF_UN_ALL;

    merge_chrom2qual(args, out);
    merge_filter(args, out);
#ifdef USE_ID_MAP
    //map from id in merged hdr to index in merged bcf1_t* out : initialize to -1
    memset(args->m_idmap.m_merged_id_2_merged_idx, -1, args->m_idmap.m_num_merged_ids*sizeof(int));
#endif
    merge_info(args, out);
#ifdef USE_ID_MAP
    //since same ID may be FORMAT also
    //map from id in merged hdr to index in merged bcf1_t* out : initialize to -1
    memset(args->m_idmap.m_merged_id_2_merged_idx, -1, args->m_idmap.m_num_merged_ids*sizeof(int));
#endif
    merge_format(args, out);

    bcf_write1(args->out_fh, args->out_hdr, out);
    if(args->m_plmedian_info.m_output_fptr)
    {
        int veclen = compute_PLmedian(args->out_hdr, out,
                &(args->m_plmedian_info.m_median_result), &(args->m_plmedian_info.m_median_result_len),
                &(args->m_plmedian_info.m_buffer), &(args->m_plmedian_info.m_buffer_len),
                &(args->m_plmedian_info.m_reorg_buffer), &(args->m_plmedian_info.m_reorg_buffer_len));
        print_PLmedian(args->m_plmedian_info.m_output_fptr, args->out_hdr, out,
                args->m_plmedian_info.m_median_result, veclen);
    }
}


void debug_buffers(FILE *fp, bcf_srs_t *files);
void debug_buffer(FILE *fp, bcf_sr_t *reader);

#define SWAP(type_t,a,b) { type_t tmp = (a); (a) = (b); (b) = tmp; }

// Clean the reader's buffer to and make it ready for the next next_line() call.
// Moves finished records (SKIP_DONE flag set) at the end of the buffer and put
// the rest to the beggining. Then shorten the buffer so that the last element
// points to the last unfinished record. There are two special cases: the last
// line of the buffer typically has a different position and must stay at the
// end; next, the first record of the buffer must be one of those already
// printed, as it will be discarded by next_line().
//
// Better explanation: For readers which participated in the last merge, reorganize buffer so that:
// a. Records [1:num_not_done] contain records with (rec->pos == pos) which are !SKIP_DONE
// b. Records [1+num_not_done:num_not_done+num_splits] contain records where the intervals are split
// c. Records [1+num_not_done+num_splits:num_diff] contain records with a different position
// d. Record [1+num_not_done+num_splits+num_diff:last_valid_record_idx]: contain records which were already processed
// Since buffer elements are fully allocated bcf1_t*, make sure you do not lose any of them (else memory leaks will occur)
// Finally update nbuffer and last_valid_record_idx
#ifdef COLLECT_STATS
void shake_buffer(args_t* args, maux_t *maux, int ir, int pos)
#else
void shake_buffer(maux_t *maux, int ir, int pos)
#endif
{
    bcf_sr_t *reader = &maux->files->readers[ir];
    bcf_hdr_t* hdr = reader->processed_header;
    maux1_t *m = maux->d[ir];

    //was not involved in merge - don't bother to check further
    if ( !reader->buffer  || reader->buffer[0]->pos != pos) return;

#define HANDLE_SPLIT_RECORD(hdr, curr_record)                                                           \
    /*new start pos is 1 greater than current END*/                                                     \
    curr_record->pos = bcf_get_end_point(curr_record) + 1;                                              \
    bcf_set_end_point(curr_record, bcf_get_split_end_point(curr_record));                               \
    /*set end point in INFO field*/                                                                     \
    bcf_set_end_point_in_info(hdr, curr_record);                                                        \
    /*set REF base*/                                                                                    \
    /*REF intervals have only 2 alleles - REF base and NON_REF tag*/                                    \
    ASSERT(curr_record->n_allele == 2);                                                                 \
    /*the REF allele in a REF interval is a single base*/                                               \
    ASSERT(strlen(curr_record->d.allele[0]) == 1);                                                      \
    curr_record->d.allele[0][0] = bcf_get_split_ref_base(curr_record);                                  \
    curr_record->d.allele[0][1] = '\0';                                                                 \
    curr_record->d.shared_dirty |= BCF1_DIRTY_ALS;                                                      \
    /*FIXME: no need to call bcf_update_alleles as the REF length is unchanged*/                        \
    /*bcf_update_alleles unnecessarily invokes free and kputc adding an overhead*/                      \
    /*bcf_update_alleles(hdr, curr_record, (const char**)curr_record->d.allele,*/                       \
    /*curr_record->n_allele);*/                                                                         \
    curr_record->m_is_split_record = 0;
    
    int i = 0;
    //Most common case, 2 entries in buffer, 1 merged in the completed round, 1 next pos
    //Common case code is as simple as possible for speed (I feel the need, the need for speed)
    if(reader->nbuffer == 1 && reader->last_valid_record_idx == 1 && (reader->buffer[1]->pos != pos)
            && (m[0].skip&SKIP_DONE))
    {
#ifdef COLLECT_STATS
        update_stat(args, NUM_FAST_SHAKE_BUFFER, 1);
#endif
        bcf1_t* tmp = 0;
        if(reader->buffer[0]->m_is_split_record)        //should re-schedule this interval
        {
            HANDLE_SPLIT_RECORD(hdr, reader->buffer[0]);
            ASSERT(reader->mbuffer >= 3);
            reader->last_valid_record_idx = reader->nbuffer = 2;
            if (maux->nbuf[ir] < 3) 
            {
                maux_expand1(maux, ir);
                m = maux->d[ir];
            }
            m[1].skip = 0;
            m[2].skip = 0;
            //shift down
            tmp = reader->buffer[2];
            reader->buffer[2] = reader->buffer[1];
            reader->buffer[1] = reader->buffer[0];
            reader->buffer[0] = tmp;
        }
        else
            reader->buffer[0]->m_is_preprocessed = 0;  //this entry is done, unset preprocessed flag
        reader->buffer[0]->pos = -1;
        return;
    }
    //Now handle the general case
#ifdef COLLECT_STATS
    update_stat(args, NUM_FAST_SHAKE_BUFFER, 0);
#endif
    // FILE *fp = stdout;
    // fprintf(fp,"<going to shake> nbuf=%d\t", reader->nbuffer); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]); fprintf(fp,"\n");
    // debug_buffer(fp,reader);
    // fprintf(fp,"--\n");
    bcf1_t** duplicate_buffer = (bcf1_t**)malloc((reader->last_valid_record_idx+1)*sizeof(bcf1_t*));
    int* duplicate_skip = (int*)malloc((reader->last_valid_record_idx+1)*sizeof(int));
    int last_same_pos_idx = (reader->buffer[reader->nbuffer]->pos == pos) ? reader->nbuffer
        : reader->nbuffer-1;
    int num_not_dones = 0;
    int num_splits = 0;
    for(i=0;i<=last_same_pos_idx;++i)
    {
        duplicate_buffer[i] = reader->buffer[i];
        duplicate_skip[i] = m[i].skip;
        if(!(m[i].skip & SKIP_DONE))
            ++num_not_dones;
        else
            if(reader->buffer[i]->m_is_split_record)
                ++num_splits;
            else
                reader->buffer[i]->m_is_preprocessed = 0;       //done and not split, unset preprocessed flag
    }
    int num_dones = (last_same_pos_idx+1)-(num_splits+num_not_dones);
    int num_diff_pos = (reader->last_valid_record_idx - last_same_pos_idx);
    for(;i<=reader->last_valid_record_idx;++i)
    {
        duplicate_skip[i] = m[i].skip;
        duplicate_buffer[i] = reader->buffer[i];
    }
    //will start writing at buffer[1] and not 0
    int not_done_idx = 1;
    int splits_idx = not_done_idx + num_not_dones;
    int diff_pos_idx = splits_idx + num_splits;
    int dones_idx = diff_pos_idx + num_diff_pos;
    int final_idx = dones_idx + num_dones - 1;
    ASSERT(final_idx == reader->last_valid_record_idx+1);
    if (final_idx >= maux->nbuf[ir] ) 
    {
        int tmp = reader->last_valid_record_idx;
        reader->last_valid_record_idx = final_idx;
        maux_expand1(maux, ir);
        reader->last_valid_record_idx = tmp;
        m = maux->d[ir];
    }
    ASSERT(reader->mbuffer > final_idx);
    //Final location will be overwritten - move to idx 0
    reader->buffer[0] = reader->buffer[final_idx];
    for(i=0;i<=last_same_pos_idx;++i)
    {
        bcf1_t* curr_record = duplicate_buffer[i];
        if(!(duplicate_skip[i]&SKIP_DONE))
        {
            m[not_done_idx].skip = duplicate_skip[i];
            reader->buffer[not_done_idx++] = curr_record;
        }
        else
        {
            if(curr_record->m_is_split_record)
            {
                HANDLE_SPLIT_RECORD(hdr, curr_record);         
                m[splits_idx].skip = 0; //should retry this record again
                reader->buffer[splits_idx++] = curr_record;
            }
            else
            {
                m[dones_idx].skip = duplicate_skip[i]; 
                reader->buffer[dones_idx++] = curr_record;
            }
        }
    }
    for(;i<=reader->last_valid_record_idx;++i)
    {
        m[diff_pos_idx].skip = duplicate_skip[i];
        reader->buffer[diff_pos_idx++] = duplicate_buffer[i];
    }
    reader->last_valid_record_idx = diff_pos_idx - 1; //ignore the dones
    if(num_not_dones > 0)
        if(num_splits + num_diff_pos == 0)      //no other valid !DONE records, point nbuffer to last !DONE record
            reader->nbuffer = not_done_idx - 1;
        else
            reader->nbuffer = not_done_idx;     //point to first record with diff pos
    else
        if(num_diff_pos == 0)                   //only splits exist, point to last split rec
            reader->nbuffer = splits_idx - 1;
        else
            reader->nbuffer = splits_idx;       //others apart from splits exist, point to first rec with diff pos
    // debug_buffer(fp,reader);
    // fprintf(fp,"<shaken>\t"); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]);
    // fprintf(fp,"\n\n");

    // set position of finished buffer[0] line to -1, otherwise swapping may
    // bring it back after next_line()
    reader->buffer[0]->pos = -1;
    free(duplicate_buffer);
    free(duplicate_skip);
}

void debug_maux(args_t *args, int pos, int var_type)
{
    bcf_srs_t *files = args->files;
    maux_t *maux = args->maux;
    int j,k,l;

    fprintf(stderr,"Alleles to merge at %d\n", pos+1);
    for (j=0; j<files->nreaders; j++)
    {
        bcf_sr_t *reader = &files->readers[j];
        fprintf(stderr," reader %d: ", j);
        for (k=0; k<=reader->nbuffer; k++)
        {
            if ( maux->d[j][k].skip==SKIP_DONE ) continue;
            bcf1_t *line = reader->buffer[k];
            if ( line->pos!=pos ) continue;
            fprintf(stderr,"\t");
            if ( maux->d[j][k].skip ) fprintf(stderr,"[");  // this record will not be merged in this round
            for (l=0; l<line->n_allele; l++)
                fprintf(stderr,"%s%s", l==0?"":",", line->d.allele[l]);
            if ( maux->d[j][k].skip ) fprintf(stderr,"]");
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr," counts: ");
    for (j=0; j<maux->nals; j++) fprintf(stderr,"%s   %dx %s", j==0?"":",",maux->cnt[j], maux->als[j]); fprintf(stderr,"\n");
    for (j=0; j<files->nreaders; j++)
    {
        bcf_sr_t *reader = &files->readers[j];
        fprintf(stderr," out %d: ", j);
        for (k=0; k<=reader->nbuffer; k++)
        {
            if ( maux->d[j][k].skip==SKIP_DONE ) continue;
            bcf1_t *line = reader->buffer[k];
            if ( line->pos!=pos ) continue;
            if ( maux->d[j][k].skip ) continue;
            fprintf(stderr,"\t");
            for (l=0; l<line->n_allele; l++)
                fprintf(stderr,"%s%s", l==0?"":",", maux->als[maux->d[j][k].map[l]]);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

//Allocates new memory
char* modify_info_string(const char* orig, char* dst, int* dst_allocated)
{
    int newlength = strlen(orig) + strlen(g_info2format_suffix)+1;
    int do_copy = 1;
    if(orig == dst)
        do_copy = 0;
    int tmp = 0;
    int* target_ptr = (dst_allocated == 0 ? &tmp : dst_allocated);
    REALLOC_IF_NEEDED(dst, *target_ptr, newlength, char);
    if(do_copy)
        strcpy(dst, orig);
    strcat(dst, g_info2format_suffix);
    return dst;
}

int get_element_size(int bcf_ht_type)
{
    switch(bcf_ht_type)
    {
        case BCF_HT_STR:
        case BCF_HT_FLAG: 
            return 1;
            break;
        case BCF_HT_INT:
            return sizeof(int);
            break;
        case BCF_HT_REAL:
            return sizeof(float);
            break;
        default:
            fprintf(stderr,"Unknown BCF_HT_TYPE : %d\n",bcf_ht_type);
    }
    assert(0);
    return -1;  //should never come here
}

void preprocess_line(args_t* args, bcf_sr_t* reader, bcf1_t* line, int file_idx)
{
    //Translate to destination header - incorrect in this case
    /*bcf_translate(reader->processed_header, reader->header, line);*/
    //Drop or move/copy INFO fields to destination's FORMAT fields
    int i = 0;
    int j = 0;
    bcf_hdr_t* src_hdr = reader->header;
    bcf_hdr_t* dst_hdr = reader->processed_header;
    //Temporary space
    args->maux->ntmp_arr = 10;
    args->maux->tmp_arr = realloc(args->maux->tmp_arr,10*sizeof(int));
    char* format_array = (char*)malloc(1024);  //Q:why 1024 A:why not? will be realloc-ed if necessary
    int format_array_allocated = 1024;
    //Array to quickly lookup actions from id
    int* id2action = g_merge_config.file2id2action[file_idx];
    //Number of samples in file
    int nsamples = bcf_hdr_nsamples(src_hdr);
    //Deal with FORMAT fields first as INFO fields may be moved/copied to FORMAT later
    for(i=0;i<line->n_fmt;++i)
    {
        bcf_fmt_t* curr_fmt = &(line->d.fmt[i]);
        int dict_id = curr_fmt->id;
        int type = bcf_hdr_id2type(src_hdr, BCF_HL_FMT, dict_id);
        //Can't declare within switch
        int delete_info = 0;
        switch(id2action[dict_id])
        {
            case MERGE_KEEP:
                break;
            case MERGE_DROP:
                //remove INFO field
                delete_info = bcf_update_format_struct(src_hdr, line, 0, 0, 0, type, curr_fmt, dict_id);
                ASSERT(delete_info >= 0);
                break;
            default:
                assert(0 && "Unknown action for format field\n");
        }
    }
    //handle INFO fields
    for(i=0;i<line->n_info;++i)
    {
        bcf_info_t* curr_info = &(line->d.info[i]);
        int dict_id = curr_info->key;
        int type = bcf_hdr_id2type(src_hdr, BCF_HL_INFO, dict_id);
        //Can't declare within switch
        int num_values_written = 0;
        int element_size = 0;
        int* flag_array = 0;
        int new_dict_id = -1;
        int idx = 0;
        int bytes_per_sample = 0;
        int delete_info = 0;
        const char* processed_info_string = 0;
        switch(id2action[dict_id])
        {
            case MERGE_KEEP:
                break;
            case MERGE_DROP:
                //remove INFO field
                delete_info = bcf_update_info_struct(src_hdr, line, 0, 0, 0, type, curr_info, dict_id);
                ASSERT(delete_info >= 0);
                break;
            case MERGE_COPY_TO_FORMAT:
            case MERGE_MOVE_TO_FORMAT:
                num_values_written = bcf_unpack_info_values(src_hdr, line, curr_info,
                        &(args->maux->tmp_arr), &(args->maux->ntmp_arr), type);
                element_size = get_element_size(type);
                //No nasty error expected
                ASSERT(num_values_written != -1 && num_values_written != -2);
                //INFO key found - move to FORMAT
                if(num_values_written > 0)
                {
                    //for flag get_info_values does not write the actual value. why? who knows
                    //Format fields have no FLAG type, convert to INT and write
                    if(type == BCF_HT_FLAG) 
                    {
                        REALLOC_IF_NEEDED(args->maux->tmp_arr, args->maux->ntmp_arr, num_values_written, int);
                        flag_array = (int*)(args->maux->tmp_arr);
                        for(j=0;j<num_values_written;++j)
                            flag_array[j] = 1;
                        element_size = sizeof(int);       //using int now
                    }
                    //Get new field name - name could be changed when moved to FORMAT as the original file could have a field with the
                    //same name
                    new_dict_id = args->m_idmap.m_readers_map[file_idx].m_original_id_2_preprocessed_copy_id[dict_id];
                    ASSERT(new_dict_id >= 0 && new_dict_id < dst_hdr->n[BCF_DT_ID]);
                    processed_info_string = dst_hdr->id[BCF_DT_ID][new_dict_id].key;
                    ASSERT(processed_info_string);
                    //either same as original info name or with _INFO suffix
                    ASSERT(strcmp(src_hdr->id[BCF_DT_ID][dict_id].key, processed_info_string) <= 0);
                    //INFO fields are common across samples, to move to FORMAT must increase the array
                    //size so that there is a value or list of values for every sample (identical info)
                    bytes_per_sample = num_values_written*element_size;
                    REALLOC_IF_NEEDED(format_array, format_array_allocated, nsamples*bytes_per_sample, char);
                    for(j=0,idx=0;j<nsamples;++j,idx+=bytes_per_sample)
                        memcpy(format_array+idx, (char*)(args->maux->tmp_arr), bytes_per_sample);
                    if(id2action[dict_id] == MERGE_MOVE_TO_FORMAT)
                    {
                        //remove INFO field
                        delete_info = bcf_update_info_struct(src_hdr, line, 0, 0, 0, type, curr_info, dict_id);
                        ASSERT(delete_info >= 0);
                    }
                    //format fields have no flag type
                    bcf_update_format(dst_hdr, line, processed_info_string, format_array, num_values_written*nsamples, type == BCF_HT_FLAG ? BCF_HT_INT : type);
#ifdef DEBUG
                    vcf_format(dst_hdr,line,&g_debug_string);
                    fprintf(g_vcf_debug_fptr,"DEBUG:preprocessed: %s",g_debug_string.s);
                    fflush(g_vcf_debug_fptr);
                    g_debug_string.l = 0;
#endif
                }
                break;
            default:
                assert(0 && "Unknown action for INFO field\n");
        }
    }
    free(format_array);
    //set preprocessed flag
    line->m_is_preprocessed = 1;
}

char get_reference_base_at_position(args_t* args, const char* seq_name, int pos)
{
    int length = 0;
    ASSERT(seq_name);
    //See if pos is within the last buffer read
    if(strcmp(args->reference_last_seq_read, seq_name) == 0 && args->reference_last_read_pos <= pos)
    {
        int offset = (pos -  args->reference_last_read_pos);
        if(offset < args->reference_num_bases_read)
            return args->reference_buffer[offset];
    }
    if(args->reference_buffer)
        free(args->reference_buffer);
    args->reference_buffer = faidx_fetch_seq(args->reference_faidx, seq_name, pos, pos+4096, &length);
    ASSERT(length > 0 && args->reference_buffer);
    strcpy(args->reference_last_seq_read, seq_name);
    args->reference_last_read_pos = pos;
    args->reference_num_bases_read = length;
    return args->reference_buffer[0];
}

//Say we are merging 3 samples S1, S2, S3
//The current records for both S1 and S2 start at position 100, end at 400 and 1000 respectively
//The next record for S3 starts at position 300
//The end point for the current interval to be merged is min (400, 1000, 300-1)
//find_end_position computes and returns this end point
//  @args : Args structure
//  @chrom_idx: Chromosome at new start of new split interval (bcf1_t.rid)
//  @ref_id: Pointer to a char* that corresponds to ID field in VCF. Filled with ID value at the start of new the split interval. Memory is
//  allocated by this function and should be freed by user
//  @ref_base: pointer to a single char, will contain the reference base at the start of new interval
int find_end_position(args_t* args, int* chrom_idx, char** ref_id, char* ref_base)
{
    bcf_srs_t* files = args->files;
    int i = 0;
    int j = 0;
    int end_point = INT32_MAX;
    int max_end_point = -1;
    //File idx not involved in current round of merge whose start pos is min
    int min_start_file_idx = -1;
    int has_line_file_idx = -1;
    int curr_end_point = 0;
    for(i=0;i<files->nreaders;++i)
    {
        bcf_sr_t* curr_reader = &(files->readers[i]);
        if(bcf_sr_has_line(files, i))
        {
            int curr_pos = curr_reader->buffer[0]->pos;
            //same pos check is necessary because of weird nature of nbuffer
            //for EOF records, buffer[nbuffer] record has same pos as others
            //in other cases, buffer[nbuffer] record has different pos as others
            int last_same_pos_idx = (curr_reader->buffer[curr_reader->nbuffer]->pos == curr_pos) ? curr_reader->nbuffer
                : curr_reader->nbuffer - 1;
            has_line_file_idx = i;
            //Multiple VCF records could have same position in a given file
            for(j=0;j<=last_same_pos_idx;++j)
            {
                bcf1_t* curr_record = curr_reader->buffer[j];
                ASSERT(curr_record && curr_record->pos >= 0 && curr_record->pos == curr_pos);
                curr_end_point = bcf_get_end_point(curr_record);
                ASSERT(curr_end_point >= 0);
                if(end_point > curr_end_point)
                {
                    end_point = curr_end_point;
                    min_start_file_idx = -1;
                }
                max_end_point = (curr_end_point > max_end_point) ? curr_end_point : max_end_point;
            }
        }
        else    //does not have the line, meaning min position record begins at greater start position
        {
            if(curr_reader->nbuffer)
            {
                bcf1_t* curr_record = curr_reader->buffer[1];   //should be min position record
                ASSERT(curr_record && curr_record->pos != -1);
                //start of this record is <= current value of end
                if(end_point >= curr_record->pos-1)
                {
                    end_point = curr_record->pos-1;
                    min_start_file_idx = i;
                }
            }
        }
    }
    ASSERT(end_point != INT32_MAX);
    //Determine base at end_point+1, i.e., the start of new interval that is created
    //If the start of new interval is the start pos of one of the records in the buffer,
    //obtain reference base etc from that record
    if(min_start_file_idx >= 0)
    {
        //Not involved in this round of merging
        ASSERT(!(bcf_sr_has_line(files, min_start_file_idx)));
        ASSERT(files->readers[min_start_file_idx].nbuffer);
        //Guaranteed to have have buffer idx 1 - obtain reference allele
        bcf1_t* rec = files->readers[min_start_file_idx].buffer[1];
        //REF + 1 variant
        ASSERT(rec->n_allele >= 2);
        //FIXME: get chromosome id and ref_id for split interval
#if 0
        //Chromosome id
        (*chrom_idx) = rec->rid;
        //ID field
        (*ref_id) = strdup(rec->d.id);
        ASSERT(*ref_id);
#endif
        //First base of reference allele - REF
        (*ref_base) = rec->d.allele[0][0];
    }
    else        //else obtain base from reference fasta file
    {
        //At least one split interval will be created, obtain REF fields from 
        if(max_end_point > end_point)
        {
            ASSERT(has_line_file_idx >= 0);
            bcf_sr_t* reader = &(files->readers[has_line_file_idx]);
            bcf1_t* rec = reader->buffer[0];
            int rid = rec->rid;
            //FIXME: get chromosome id and ref_id for split interval
#if 0
            //Chromosome id
            (*chrom_idx) = rid;
            //ID field - set to empty (NULL) since we don't know the ID at new interval's position
            (*ref_id) = 0;
#endif
            //First base of reference allele - REF, obtain from reference file
            const char* contig_name = reader->processed_header->id[BCF_DT_CTG][rid].key;
            (*ref_base) = get_reference_base_at_position(args, contig_name, end_point+1);
        }
    }
    return end_point;
}

void store_split_info(bcf_hdr_t* curr_header, bcf1_t* curr_record, int new_end_point,
        int chrom_idx, const char* ref_id, char ref_base)
{
    curr_record->m_is_split_record = 1;
    //old end point will be end point of split interval
    bcf_set_split_end_point(curr_record, bcf_get_end_point(curr_record));
    //Set REF base at start of new split interval
    bcf_set_split_ref_base(curr_record, ref_base);
    //update end point
    bcf_set_end_point(curr_record, new_end_point);
    //FIXME: update ref_id, chrom_idx
#if 0
    bcf_set_split_reference_base(curr_record, ref_base);
    //Update chrom idx
    duplicate->rid = chrom_idx;
    //Update ref_base
    
    return duplicate;
#endif
}

//Say we are merging 3 samples S1, S2, S3
//The current records for both S1 and S2 start at position 100, end at 400 and 1000 respectively
//The next record for S3 starts at position 300
//The end point for the current interval to be merged is min (400, 1000, 300-1)
//align_to_end aligns end points of all records to merge in this round to end_point
//For our example, the intervals of S1 and S2 to be merged now become (100, 299)
//Additional bcf records are created for intervals (300-400) for S1 and (300-1000) for S2
//These additional records are stored at the end of the buffer - since mbuffer >= nbuffer+1
//S3 is not touched
void align_to_end(args_t* args, int new_end_point, int chrom_idx, const char* ref_id, char ref_base)
{
    int i = 0;
    int j = 0;
    bcf_srs_t* files = args->files;
    int end_point = -1;
    for(i=0;i<files->nreaders;++i)
    {
        bcf_sr_t* curr_reader = &(files->readers[i]);
        bcf_hdr_t* curr_header = curr_reader->processed_header;
        //Only files with the min start position are considered
        if(bcf_sr_has_line(files, i))
        {
            int curr_pos = curr_reader->buffer[0]->pos;
            //same pos check is necessary because of weird nature of nbuffer
            //for EOF records, buffer[nbuffer] record has same pos as others
            //in other cases, buffer[nbuffer] record has different pos as others
            int last_same_pos_idx = (curr_reader->buffer[curr_reader->nbuffer]->pos == curr_pos) ? curr_reader->nbuffer
                : curr_reader->nbuffer - 1;
            int num_split_records_added = 0;
            for(j=0;j<=last_same_pos_idx;++j)
            {
                bcf1_t* curr_record = curr_reader->buffer[j];
                ASSERT(curr_record && curr_record->pos >= 0);
                end_point = bcf_get_end_point(curr_record);
                ASSERT(end_point >= new_end_point);
                if(end_point > new_end_point)      //need to split record
                {
                    //store information about the split
                    store_split_info(curr_header, curr_record, new_end_point, chrom_idx, ref_id, ref_base);
                    ++num_split_records_added;
                }
            }
        }
    }
}

// Determine which line should be merged from which reader: go through all
// readers and all buffered lines, expand REF,ALT and try to match lines with
// the same ALTs. A step towards output independent on input ordering of the
// lines.
void merge_buffer(args_t *args)
{
    bcf_srs_t *files = args->files;
    char *id = NULL;
    maux_t *maux = args->maux;
    //Same buffer could produce lots of lines for gVCFs
    //Invariant at the beginning of every iteration:
    // buffer[0] contains valid record with min. position for all samples which have a record at that position, invalid = -1
    // reader->nbuffer - Confusing usage
    // The best way to remember is that it is the idx in reader->buffer at which the last record read from file is stored
    // In the normal case, the pos of buffer[nbuffer] is different from all the previous records is located.
    // The exception is when EOF is reached, at that point buffer[nbuffer] has same pos as previous records
    //    case 1: reader has a valid record at idx 0, implying reader contains a min-position record to merge now
    //       [0:nbuffer-1]: valid records with same pos 
    //       [nbuffer]: valid record with different pos value
    //       The number of valid entries is nbuffer+1
    //    case 2: reader has no valid entry at idx 0, reader contains no record to be merged now
    //       0: invalid BCF record
    //       [1:nbuffer-1]: valid records with same pos
    //       [nbuffer]: valid record with different post value
    //       The number of valid entries is nbuffer

    int i, pos = -1, var_type = 0;
    maux_reset(maux);

    //Drop or move fields if necessary
    if(g_preprocess_vcfs)
    {
        for (i=0; i<files->nreaders; i++)
            if ( bcf_sr_has_line(files,i) )
            {
                bcf_sr_t* curr_reader = &(files->readers[i]);
                int j = 0;
                int last_same_pos_idx = (curr_reader->buffer[0]->pos == 
                        curr_reader->buffer[curr_reader->nbuffer]->pos)
                    ? curr_reader->nbuffer : curr_reader->nbuffer-1;
                for(j=0;j<=last_same_pos_idx;++j)
                {
                    bcf1_t* line = curr_reader->buffer[j];
                    bcf_unpack(line, BCF_UN_ALL);
                    if(!line->m_is_preprocessed)
                        preprocess_line(args, curr_reader, line, i);
                }
            }
    }

    int end_point = -1;
    //Find the 'end' point where interval has to 'break' in gVCFs
    //Create intervals with the end position
    //TODO: compute end points once before the while loop using sorting
    if(g_is_input_gvcf)
    {
        int chrom_idx = -1;
        char* ref_id = 0;
        char ref_base = 0;
        //Do unpack and don't bother to check later
        for (i=0; i<files->nreaders; i++)
            if ( bcf_sr_has_line(files,i) )
            {
                int j = 0;
                int last_same_pos_idx = (files->readers[i].buffer[0]->pos == 
                        files->readers[i].buffer[files->readers[i].nbuffer]->pos)
                    ? files->readers[i].nbuffer : files->readers[i].nbuffer-1;
                for(j=0;j<=last_same_pos_idx;++j)
                {
                    bcf_unpack(files->readers[i].buffer[j], BCF_UN_INFO);
                    bcf_set_end_point_from_info(files->readers[i].processed_header, files->readers[i].buffer[j]);
                }
            }
        end_point = find_end_position(args, &chrom_idx, &ref_id, &ref_base);
        align_to_end(args, end_point, chrom_idx, ref_id, ref_base);
        if(ref_id)
            free(ref_id);
        for (i=0; i<files->nreaders; i++)
            if ( bcf_sr_has_line(files,i) )
            {
                int j = 0;
                int last_same_pos_idx = (files->readers[i].buffer[0]->pos == 
                        files->readers[i].buffer[files->readers[i].nbuffer]->pos)
                    ? files->readers[i].nbuffer : files->readers[i].nbuffer-1;
                for(j=0;j<=last_same_pos_idx;++j)
                    bcf_set_end_point_in_info(files->readers[i].processed_header, files->readers[i].buffer[j]);
            }
    }

    // set the current position
    for (i=0; i<files->nreaders; i++)
    {
        if ( bcf_sr_has_line(files,i) )
        {
            bcf1_t *line = bcf_sr_get_line(files,i);
            pos = line->pos;
            var_type = bcf_get_variant_types(line);
            id = line->d.id;
            break;
        }
    }
    if(g_measure_iterator_timing_only)
    {
        for (i=0; i<files->nreaders; i++)
        {
            bcf_sr_t *reader = &files->readers[i];
            if ( !reader->buffer ) continue;
            int j = 0;
            for (j=0; j<=reader->nbuffer; j++)  
            {
                bcf1_t *line = reader->buffer[j];
                if(line->pos == pos)
                    maux->d[i][j].skip |= SKIP_DONE;
                else
                    break;
            }
        }
    }
    else
    {
        int contains_non_ref = 0;
        // In this loop we select from each reader compatible candidate lines.
        // (i.e. SNPs or indels). Go through all files and all lines at this
        // position and normalize relevant alleles.
        // REF-only sites may be associated with both SNPs and indels.
        for (i=0; i<files->nreaders; i++)
        {
            bcf_sr_t *reader = &files->readers[i];
            if ( !reader->buffer ) continue;
            int j, k;
            for (j=0; j<=reader->nbuffer; j++)  
            {
                bcf1_t *line = reader->buffer[j];
                int line_type = bcf_get_variant_types(line);

                // select relevant lines
                maux->d[i][j].skip = SKIP_DIFF;
                if ( pos!=line->pos ) 
                {
                    if ( j==0 ) maux->d[i][j].skip |= SKIP_DONE; // left from previous run, force to ignore
                    continue; 
                }
                if ( args->merge_by_id )
                {
                    if ( strcmp(id,line->d.id) ) continue;
                }
                else
                {
                    if ( args->collapse==COLLAPSE_NONE && maux->nals )
                    {
                        // All alleles of the tested record must be present in the
                        // selected maux record plus variant types must be the same
                        if ( var_type!=line->d.var_type ) continue;
                        if ( vcmp_set_ref(args->vcmp,maux->als[0],line->d.allele[0]) < 0 ) continue;   // refs not compatible
                        for (k=1; k<line->n_allele; k++)
                        {
                            if ( vcmp_find_allele(args->vcmp,maux->als+1,maux->nals-1,line->d.allele[k])>=0 ) break;
                        }
                        if ( k==line->n_allele ) continue;  // no matching allele
                    }
                    if ( var_type&VCF_SNP && !(line_type&VCF_SNP) && !(args->collapse&COLLAPSE_ANY) && line_type!=VCF_REF ) continue;
                    if ( var_type&VCF_INDEL && !(line_type&VCF_INDEL) && !(args->collapse&COLLAPSE_ANY) && line_type!=VCF_REF ) continue;
                }
                maux->d[i][j].skip = 0;

                //realloc
                hts_expand(int, line->n_allele, maux->d[i][j].mmap, maux->d[i][j].map);
                if ( !maux->nals )    // first record, copy the alleles to the output
                {
                    maux->nals = 0;
                    hts_expand0(char*, line->n_allele, maux->mals, maux->als);
                    hts_expand0(int, line->n_allele, maux->ncnt, maux->cnt);
                    for (k=0; k<line->n_allele; k++)
                    {
                        //NON_REFS will be added at the end
                        if(line->d.var[k].type & VCF_NON_REF)
                        {
                            contains_non_ref = 1;
                            continue;
                        }
                        maux->als[maux->nals] = strdup(line->d.allele[k]);
                        maux->d[i][j].map[k] = maux->nals;
                        maux->cnt[maux->nals] = 1;
                        ++(maux->nals);
                    }
                    pos = line->pos;
                    continue;
                }

                // normalize alleles
                maux->als = merge_alleles(line->d.allele, line->n_allele, maux->d[i][j].map, maux->als, &maux->nals, &maux->mals,
                        line->d.var, &contains_non_ref, line->m_is_split_record);
                if ( !maux->als ) error("Failed to merge alleles at %s:%d\n",bcf_seqname(args->out_hdr,line),line->pos+1);
                hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
                for (k=1; k<line->n_allele; k++)
                {
                    if(line->d.var[k].type & VCF_NON_REF)
                        continue;
                    maux->cnt[ maux->d[i][j].map[k] ]++;    // how many times an allele appears in the files
                }
                maux->cnt[0]++;
            }
        }
        //Add NON_REF variant at the end
        if(contains_non_ref)
        {
            ++(maux->nals);
            hts_expand0(char*, maux->nals, maux->mals, maux->als);
            hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
            int non_ref_idx = maux->nals - 1;
            maux->als[non_ref_idx] = strdup("<NON_REF>"); 
            for (i=0; i<files->nreaders; i++)
            {
                bcf_sr_t *reader = &files->readers[i];
                if ( !reader->buffer ) continue;
                int j, k;
                for (j=0; j<=reader->nbuffer; j++)
                {
                    bcf1_t *line = reader->buffer[j];
                    if ( pos!=line->pos || maux->d[i][j].skip ) 
                        continue;
                    //don't bother with reference allele
                    for (k=1; k<line->n_allele; k++)
                    {
                        if(line->d.var[k].type & VCF_NON_REF)
                        {
                            maux->d[i][j].map[k] = non_ref_idx;
                            ++(maux->cnt[non_ref_idx]);
                        }
                    }
                }
            }
        }
        // debug_maux(args, pos, var_type);

        // Select records that have the same alleles; the input ordering of indels
        // must not matter. Multiple VCF lines can be emitted from this loop.
        // We expect only very few alleles and not many records with the same
        // position in the buffers, therefore the nested loops should not slow us
        // much.
        while (1)
        {
            // take the most frequent allele present in multiple files
            int icnt = 0;
            for (i=1; i<maux->nals; i++) 
                if ( maux->cnt[i] > maux->cnt[icnt] ) icnt = i;
            if ( maux->cnt[icnt]<0 ) break;

            int nmask = 0;
            for (i=0; i<files->nreaders; i++)
            {
                maux->has_line[i] = 0;

                bcf_sr_t *reader = &files->readers[i];
                if ( !reader->buffer ) continue;

                // find lines with the same allele
                int j;
                for (j=0; j<=reader->nbuffer; j++)
                {
                    if ( maux->d[i][j].skip ) continue;
                    int k;
                    for (k=0; k<reader->buffer[j]->n_allele; k++)
                        if ( icnt==maux->d[i][j].map[k] ) break;
                    if ( k<reader->buffer[j]->n_allele ) break;
                }
                if ( j>reader->nbuffer )
                {
                    // no matching allele found in this file 
                    if ( args->collapse==COLLAPSE_NONE ) continue;

                    for (j=0; j<=reader->nbuffer; j++)
                    {
                        if ( maux->d[i][j].skip ) continue;
                        if ( args->collapse&COLLAPSE_ANY ) break;
                        int line_type = bcf_get_variant_types(reader->buffer[j]);
                        if ( var_type&VCF_SNP && line_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                        if ( var_type&VCF_INDEL && line_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                        if ( line_type==VCF_REF )
                        {
                            if ( var_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                            if ( var_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                        }
                        else if ( var_type==VCF_REF )
                        {
                            if ( line_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                            if ( line_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                        }
                    }
                }
                if ( j<=reader->nbuffer ) 
                {
                    // found a suitable line for merging, place it at the beggining
                    if ( j>0 ) 
                    {
                        SWAP(bcf1_t*, reader->buffer[0], reader->buffer[j]); 
                        SWAP(maux1_t, maux->d[i][0], maux->d[i][j]); 
                    }
                    // mark as finished so that it's ignored next time
                    maux->d[i][0].skip |= SKIP_DONE;
                    maux->has_line[i] = 1;
                    nmask++;
                }
            }
            if ( !nmask ) break;    // done, no more lines suitable for merging found 
#ifdef COLLECT_STATS
            update_stat(args, NUM_ACTIVE_SAMPLES_PER_LINE, nmask);
#endif
            merge_line(args);       // merge and output the line
            maux->cnt[icnt] = -1;   // do not pick this allele again, mark it as finished
        }

        // clean the alleles
        for (i=0; i<maux->nals; i++)
        {
            free(maux->als[i]);
            maux->als[i] = 0;
        }
        maux->nals = 0;
    }

    // get the buffers ready for the next next_line() call
    for (i=0; i<files->nreaders; i++)
#ifdef COLLECT_STATS
        shake_buffer(args, maux, i, pos);
#else
        shake_buffer(maux, i, pos);
#endif
}

void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd)
{
    kstring_t str = {0,0,0};
    ksprintf(&str,"##%sVersion=%s+htslib-%s\n", cmd, bcftools_version(), hts_version());
    bcf_hdr_append(hdr,str.s);

    str.l = 0;
    ksprintf(&str,"##%sCommand=%s", cmd, argv[0]);
    int i;
    for (i=1; i<argc; i++)
    {
        if ( strchr(argv[i],' ') )
            ksprintf(&str, " '%s'", argv[i]);
        else
            ksprintf(&str, " %s", argv[i]);
    }
    kputc('\n', &str);
    bcf_hdr_append(hdr,str.s);
    free(str.s);

    bcf_hdr_sync(hdr);
}

typedef struct
{
    int type; //BCF_HL_*
    int id;
    const char* key;
} erase_hrec_struct;

typedef struct
{
    int m_original_id;
    char* m_new_key;
    bcf_hrec_t* m_new_hrec;
} add_hrec_struct;

//Drop unnecessary fields
//Move/copy INFO fields to FORMAT
void preprocess_header(args_t* args, int file_idx)
{
    int i = 0;
    int j = 0;
    bcf_hdr_t* hdr = args->files->readers[file_idx].processed_header;

    //Array to store hrecs which will be erased from hdr
    erase_hrec_struct* erase_array = (erase_hrec_struct*)malloc(hdr->n[BCF_DT_ID]*sizeof(erase_hrec_struct));
    int erase_idx = 0;
    //Array to store hrecs which will be added at the end
    add_hrec_struct* add_array = (add_hrec_struct*)malloc(hdr->n[BCF_DT_ID]*sizeof(add_hrec_struct));
    int add_idx = 0;

    //Used in preprocess_line to quickly lookup action from id (integer id)
    ASSERT(file_idx < g_merge_config.num_files);
    g_merge_config.file2id2action[file_idx] = (int*)malloc(hdr->n[BCF_DT_ID]*sizeof(int));
    //Drop fields which need to be dropped
    //Move INFO fields to FORMAT
    for(i=0;i<hdr->n[BCF_DT_ID];++i)        //iterate over ID key-value pairs
    {
        //is this an INFO or FORMAT field?
        int is_info = bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, i);
        int is_format = bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, i);
        g_merge_config.file2id2action[file_idx][i] = g_merge_config.default_action;
        //Map from id to id - same as duplicate
        args->m_idmap.m_readers_map[file_idx].m_original_id_2_preprocessed_id[i] = i;
        if(is_info || is_format)
        {
            bcf_idpair_t* curr_pair = &(hdr->id[BCF_DT_ID][i]);
            const char* key_string = curr_pair->key;
            int action_type = get_merge_action(key_string, is_info, is_format);
            //Used in preprocess_line to quickly lookup action from id (integer id)
            g_merge_config.file2id2action[file_idx][i] = action_type;
            //END in INFO has special significance
            if(is_info && strcmp(key_string,"END") == 0)        
            {
                g_merge_config.file2id2action[file_idx][i] = MERGE_KEEP;  //keep END in INFO always
                continue;
            }
            int type = is_info ? bcf_hdr_id2type(hdr, BCF_HL_INFO, i) : bcf_hdr_id2type(hdr, BCF_HL_FMT, i);
            // Is the action drop or move?
            if(do_detach(action_type))
            {
                //Mark for removal
                ASSERT(erase_idx < hdr->n[BCF_DT_ID]);
                erase_array[erase_idx].id = i;
                erase_array[erase_idx].type = is_info ? BCF_HL_INFO : BCF_HL_FMT;
                erase_array[erase_idx++].key = key_string;
                //This field does not exist, mark invalid
                args->m_idmap.m_readers_map[file_idx].m_original_id_2_preprocessed_id[i] = -1;
            }
            if(is_format && g_do_gatk_merge && strcmp(key_string,"AD") == 0)  //For AD, Number=R - in some gVCFs the Number is set to '.'
                bcf_hdr_set_length_descriptor(hdr, BCF_HL_FMT, i, BCF_VL_R);
            //If move to format, add to FORMAT
            if(is_info && (action_type == MERGE_MOVE_TO_FORMAT || action_type == MERGE_COPY_TO_FORMAT))
            {
                //if there is a FORMAT field with the exact same name, don't bother to move/copy
                //instead just DROP if move is requested
                if(bcf_hdr_field_type_id2int(hdr,key_string, BCF_DT_ID, BCF_HL_FMT) >= 0)
                {
                    if(action_type == MERGE_MOVE_TO_FORMAT)
                        g_merge_config.file2id2action[file_idx][i] = MERGE_DROP;
                    else
                        g_merge_config.file2id2action[file_idx][i] = MERGE_KEEP;
                    continue;
                }
                //Create a copy of hrec as we will remove that hrec from the header and add new hrec
                bcf_hrec_t* orig_hrec = bcf_hdr_id2hrec(hdr, BCF_DT_ID, BCF_HL_INFO, i);
                ASSERT(orig_hrec);
                bcf_hrec_t* new_hrec = bcf_hrec_dup(orig_hrec);
                //change type to FORMAT
                new_hrec->type = BCF_HL_FMT;
                new_hrec->key = (char*)realloc(new_hrec->key, strlen("FORMAT")+1);        //for end NULL char
                strcpy(new_hrec->key, "FORMAT");
                //Find ID field in INFO hrec and change its value
                int idx = bcf_hrec_find_key(new_hrec, "ID");
                ASSERT(idx >= 0);
                //If there is a FORMAT field with the same name, add _INFO suffix
                char* new_key_string = new_hrec->vals[idx];
                /*modify_info_string(new_hrec->vals[idx], new_hrec->vals[idx], 0) : */
                new_hrec->vals[idx] = new_key_string;
                //FORMAT fields do not have FLAG type, convert to int if this field is of type FLAG
                if(type == BCF_HT_FLAG)
                {
                    //Number == 0 for all FLAG types in INFO, set to 1 for FORMAT field
                    idx = bcf_hrec_find_key(new_hrec, "Number");
                    if(idx <  0)    //No such field exists, add field first
                    {
                        bcf_hrec_add_key(new_hrec,"Number",strlen("Number")); 
                        idx = bcf_hrec_find_key(new_hrec, "Number");
                    }
                    ASSERT(idx >= 0);
                    bcf_hrec_set_val(new_hrec, idx, "1", 1, 0);
                    //Change type to Integer
                    idx = bcf_hrec_find_key(new_hrec, "Type");
                    if(idx <  0)    //No such field exists, add field first
                    {
                        bcf_hrec_add_key(new_hrec,"Type",strlen("Type")); 
                        idx = bcf_hrec_find_key(new_hrec, "Type");
                    }
                    ASSERT(idx >= 0);
                    bcf_hrec_set_val(new_hrec, idx, "Integer", strlen("Integer"), 0);
                }
                ASSERT(add_idx < hdr->n[BCF_DT_ID]);
                //Store format field to add to new hdr
                add_array[add_idx].m_original_id = i;
                add_array[add_idx].m_new_key = new_key_string;
                add_array[add_idx++].m_new_hrec = new_hrec;
            }
        }
    }
    //Erase hrecs
    for(j=0;j<erase_idx;++j)
        bcf_hdr_remove(hdr, BCF_DT_ID, erase_array[j].key);
    free(erase_array);
    //Add hrecs
    for(j=0;j<add_idx;++j)
    {
        bcf_hdr_add_hrec(hdr, add_array[j].m_new_hrec);
        int new_id = bcf_hdr_field_type_id2int(hdr, add_array[j].m_new_key, BCF_DT_ID, add_array[j].m_new_hrec->type);
        ASSERT(new_id >= 0);
        //Store mapping from original id to new id
        //The field is copied - so store in copy id array
        int original_id = add_array[j].m_original_id;
        args->m_idmap.m_readers_map[file_idx].m_original_id_2_preprocessed_copy_id[original_id] = new_id;
    }
    free(add_array);
#ifdef DEBUG
    int hdr_str_length = 0;
    char* hdr_str = bcf_hdr_fmt_text(hdr, 0, &hdr_str_length);
    fprintf(g_debug_fptr,"%s",hdr_str);
    fflush(g_debug_fptr);
    free(hdr_str);
#endif
}

void allocate_original_idmap(args_t* args)
{
    int i = 0;
    //map from ids of file i to merged hdr ids
    args->m_idmap.m_readers_map = (reader_idmap*)calloc(args->files->nreaders, sizeof(reader_idmap));
    for(i=0;i<args->files->nreaders;++i)
    {
        reader_idmap* curr_map = &(args->m_idmap.m_readers_map[i]);
        bcf_hdr_t* hdr = args->files->readers[i].header;
        curr_map->m_num_original_ids = hdr->n[BCF_DT_ID];
        curr_map->m_original_id_2_preprocessed_id = (int*)malloc(curr_map->m_num_original_ids*sizeof(int));
        curr_map->m_original_id_2_preprocessed_copy_id = (int*)malloc(curr_map->m_num_original_ids*sizeof(int));
        //set to invalid value -1
        memset(curr_map->m_original_id_2_preprocessed_id, -1, curr_map->m_num_original_ids*sizeof(int));
        memset(curr_map->m_original_id_2_preprocessed_copy_id, -1, curr_map->m_num_original_ids*sizeof(int));
    }
}

void setup_idmap(args_t* args)
{
    int i = 0;
    int j = 0;
    //map from ids of file i to merged hdr ids
    if(args->m_idmap.m_readers_map == 0)
        allocate_original_idmap(args);
    bcf_hdr_t* out_hdr = args->out_hdr;
    args->m_idmap.m_num_merged_ids = out_hdr->n[BCF_DT_ID];
    //map from merged ids to ids of file i
    args->m_idmap.m_merged_id_2_reader_idmap = (int**)malloc(args->files->nreaders*sizeof(int*));
    //map from id in merged hdr to index in merged bcf1_t* out
    args->m_idmap.m_merged_id_2_merged_idx = (int*)malloc(args->m_idmap.m_num_merged_ids*sizeof(int));
    //map from id in merged hdr to info_rule index
    args->m_idmap.m_merged_id_2_info_rule_idx = (int*)malloc(args->m_idmap.m_num_merged_ids*sizeof(int));
    memset(args->m_idmap.m_merged_id_2_info_rule_idx, -1, args->m_idmap.m_num_merged_ids*sizeof(int));
    args->m_idmap.m_DP_info_vals = (int*)malloc(args->files->nreaders*sizeof(int));
    for(i=0;i<args->files->nreaders;++i)
    {
        //map from merged ids to ids of file i
        args->m_idmap.m_merged_id_2_reader_idmap[i] = (int*)malloc(out_hdr->n[BCF_DT_ID]*sizeof(int));
        //initialize to -1
        memset(args->m_idmap.m_merged_id_2_reader_idmap[i], -1, out_hdr->n[BCF_DT_ID]*sizeof(int));
        bcf_hdr_t* curr_hdr = args->files->readers[i].processed_header;
        //map from ids of file i to merged hdr ids
        reader_idmap* curr_map = &(args->m_idmap.m_readers_map[i]); 
        curr_map->m_num_ids = curr_hdr->n[BCF_DT_ID];
        curr_map->m_id_2_merged_id = (int*)malloc(curr_hdr->n[BCF_DT_ID]*sizeof(int));
        curr_map->m_GQ_id = bcf_hdr_id2int(curr_hdr, BCF_DT_ID, "GQ");
        for(j=0;j<curr_hdr->n[BCF_DT_ID];++j)
        {
            bcf_idpair_t* curr_id = &(curr_hdr->id[BCF_DT_ID][j]);
            if(curr_id->key == 0)       //deleted
                continue;
            ASSERT(curr_id->val->id == j);
            int col_type = bcf_hdr_idinfo_exists(curr_hdr,BCF_HL_FLT, j) ? BCF_HL_FLT : 
                bcf_hdr_idinfo_exists(curr_hdr,BCF_HL_INFO, j) ? BCF_HL_INFO :
                bcf_hdr_idinfo_exists(curr_hdr,BCF_HL_FMT, j) ? BCF_HL_FMT : -1;
            ASSERT(col_type != -1);
            bcf_hrec_t* hrec = bcf_hdr_id2hrec(curr_hdr, BCF_DT_ID, col_type, j);
            int out_idx = -1;
            if(hrec)       //not deleted
            {
                const char* key = curr_id->key;
                out_idx = bcf_hdr_id2int(out_hdr, BCF_DT_ID, key);
#ifdef DEBUG
                ASSERT(out_idx >= 0 && out_idx < out_hdr->n[BCF_DT_ID]);
                ASSERT(bcf_hdr_idinfo_exists(out_hdr, col_type, out_idx));  //same coltype as input
#endif
                args->m_idmap.m_merged_id_2_reader_idmap[i][out_idx] = j;
            }
            curr_map->m_id_2_merged_id[j] = out_idx;
        }
    }
    args->m_idmap.m_merged_DP_info_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "DP");
    args->m_idmap.m_merged_DP_format_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "DP");
    args->m_idmap.m_merged_MIN_DP_format_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "MIN_DP");
    ASSERT(args->m_idmap.m_merged_DP_format_id == args->m_idmap.m_merged_DP_info_id);
    ASSERT(args->m_idmap.m_merged_DP_info_id < 0 || bcf_hdr_idinfo_exists(out_hdr, BCF_HL_INFO, args->m_idmap.m_merged_DP_info_id));
    ASSERT(args->m_idmap.m_merged_DP_info_id < 0 || bcf_hdr_idinfo_exists(out_hdr, BCF_HL_FMT, args->m_idmap.m_merged_DP_info_id));
    args->m_idmap.m_merged_END_info_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "END");
}

void destroy_idmap(args_t* args)
{
    int i = 0;
    if(args->m_idmap.m_readers_map)
    {
        for(i=0;i<args->files->nreaders;++i)
        {
            reader_idmap* ptr = &(args->m_idmap.m_readers_map[i]);
            free(ptr->m_id_2_merged_id);
            free(ptr->m_original_id_2_preprocessed_id);
            free(ptr->m_original_id_2_preprocessed_copy_id);
            free(args->m_idmap.m_merged_id_2_reader_idmap[i]);
        }
        free(args->m_idmap.m_readers_map);
        free(args->m_idmap.m_merged_id_2_reader_idmap);
        free(args->m_idmap.m_merged_id_2_merged_idx);
        free(args->m_idmap.m_merged_id_2_info_rule_idx);
        free(args->m_idmap.m_DP_info_vals);
    }
}

void merge_vcf(args_t *args)
{
    args->out_fh  = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    args->out_hdr = bcf_hdr_init("w");
    //set collapse - BUG in original code; needed for sync_reader to work correctly
    args->files->collapse = args->collapse;

    if ( args->header_fname )
    {
        if ( bcf_hdr_set(args->out_hdr,args->header_fname) ) error("Could not read/parse the header: %s\n", args->header_fname);
    }
    else
    {
        int i;
        if(g_preprocess_vcfs)
        {
            g_merge_config.file2id2action = (int**)malloc(args->files->nreaders*sizeof(int*));
            g_merge_config.num_files = args->files->nreaders;
        }
        allocate_original_idmap(args);
        for (i=0; i<args->files->nreaders; i++)
        {
            char buf[10]; snprintf(buf,10,"%d",i+1);
            bcf_hdr_t* hdr = args->files->readers[i].header;
            if(g_preprocess_vcfs)       //create duplicate header for changes
            {
                args->files->readers[i].header->ntransl = 0;    //required for bcf_translate later
                args->files->readers[i].processed_header = bcf_hdr_dup(args->files->readers[i].header);
                preprocess_header(args, i);
                hdr = args->files->readers[i].processed_header;
            }
            else                //else point to same header
                args->files->readers[i].processed_header = args->files->readers[i].header;
            bcf_hdr_merge(args->out_hdr, hdr, buf, args->force_samples);
        }
        bcf_hdr_append_version(args->out_hdr, args->argc, args->argv, "bcftools_merge");
        bcf_hdr_sync(args->out_hdr);
        setup_idmap(args);
    }
    info_rules_init(args);

    bcf_hdr_set_version(args->out_hdr, bcf_hdr_get_version(args->files->readers[0].header));
    bcf_hdr_write(args->out_fh, args->out_hdr);
    if ( args->header_only )
    {
        bcf_hdr_destroy(args->out_hdr);
        hts_close(args->out_fh);
        return;
    }

    if ( args->collapse==COLLAPSE_NONE ) args->vcmp = vcmp_init();
    args->maux = maux_init(args->files);
    args->out_line = bcf_init1();
    args->tmph = kh_init(strdict);
    int ret;
#ifdef PROFILE
    ProfilerStart("gprofile.log");
#endif
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        merge_buffer(args);
    }
#ifdef PROFILE
    ProfilerStop();
#endif
    info_rules_destroy(args);
    maux_destroy(args->maux);
    bcf_hdr_destroy(args->out_hdr);
    hts_close(args->out_fh);
    bcf_destroy1(args->out_line);
    kh_destroy(strdict, args->tmph);
    if ( args->tmps.m ) free(args->tmps.s);
    if ( args->vcmp ) vcmp_destroy(args->vcmp);
}

void initialize_reference(args_t* args, const char* reference_filename)
{
    args->reference_filename = strdup(reference_filename);
    assert(args->reference_filename && "char* storing reference file name in args is NULL");
    /*args->reference_fp = gzopen(argv[1], "r");*/
    /*assert(args->reference_fp != Z_NULL && "Failed to open reference file");*/
    /*args->reference_seq = kseq_init(args->reference_fp);*/
    /*assert(args->reference_seq && "Unable to allocate memory for reference seq");*/
    args->reference_faidx = fai_load(args->reference_filename);
    assert(args->reference_faidx);
    args->reference_last_seq_read = (char*)calloc(4096, sizeof(char));
}

void destroy_reference(args_t* args)
{
    if(args->reference_filename)
    {
        fai_destroy(args->reference_faidx);
        free(args->reference_filename);
        free(args->reference_last_seq_read);
        if(args->reference_buffer)
            free(args->reference_buffer);
    }
}

void parse_merge_config(char* filename)
{
    FILE* fptr = fopen(filename,"r");
    if(fptr == 0)
    {
        fprintf(stderr,"Unable to open merge config file %s -- exiting\n",filename);
        exit(-1);
    }
    g_merge_config.is_initialized = 1;
    g_preprocess_vcfs = 1;
    int ret = 0;
    khiter_t kitr;
    //Initialize tokens (char*) to idx - so that in future we use int instead of char*
    g_merge_config_token_2_idx = kh_init(strdict);
    kh_clear(strdict, g_merge_config_token_2_idx);
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "DEFAULT", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_DEFAULT;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "INFO", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_INFO;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "FORMAT", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_FORMAT;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "KEEP", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_KEEP;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "DROP", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_DROP;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "MOVE_TO_FORMAT", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_MOVE_TO_FORMAT;
    kitr = kh_put(strdict, g_merge_config_token_2_idx, "COPY_TO_FORMAT", &ret);
    kh_val(g_merge_config_token_2_idx, kitr) = MERGE_COPY_TO_FORMAT;

    //Initialize merge_config structure
    g_merge_config.default_action = MERGE_KEEP;
    g_merge_config.info_field2action = kh_init(strdict);
    kh_clear(strdict, g_merge_config.info_field2action);
    g_merge_config.format_field2action = kh_init(strdict);
    kh_clear(strdict, g_merge_config.format_field2action);

    char* line = 0;
    size_t linesize = 0;
    size_t length = 0;
    const char* delim = " \t";
    int line_number = 1;
    while(1)
    {
        length = getline(&line, &linesize, fptr);
        if(feof(fptr))
            break;
        if(length > 0 && line[length-1] == '\n')
        {
            --length;
            line[length] = '\0';
        }
        if(length == 0)
            continue;
        char* ptr = strtok(line, delim);
        if(ptr == 0 || ptr[0] == '#')       //comment
            continue;
        kitr = kh_get(strdict, g_merge_config_token_2_idx, ptr);
        if( kitr == kh_end(g_merge_config_token_2_idx) )
        {
            fprintf(stderr,"Unknown token %s in line %d in merge config file %s\n", ptr, line_number, filename);
            exit(-1);
        }
        int token_type = kh_val(g_merge_config_token_2_idx, kitr);
        char* field = 0;
        switch(token_type)
        {
            case MERGE_DEFAULT:
                ptr = strtok(0, delim);
                assert(ptr && "DEFAULT token needs an argument : KEEP|DROP");
                kitr = kh_get(strdict, g_merge_config_token_2_idx, ptr);
                assert(kitr != kh_end(g_merge_config_token_2_idx) && "Argument to DEFAULT must be KEEP|DROP");
                g_merge_config.default_action = kh_val(g_merge_config_token_2_idx, kitr);
                assert((g_merge_config.default_action == MERGE_KEEP || g_merge_config.default_action == MERGE_DROP)
                    && "Argument to DEFAULT must be KEEP|DROP");
                break;
            case MERGE_INFO:
            case MERGE_FORMAT:
                field = strdup(strtok(0, delim));
                ptr = strdup(strtok(0, delim)); //action
                assert(field && ptr && "INFO/FORMAT token needs 2 arguments : FIELD_NAME KEEP|DROP etc.");
                kitr = kh_get(strdict, g_merge_config_token_2_idx, ptr);
                assert(kitr != kh_end(g_merge_config_token_2_idx) && "Unknown action for INFO/FORMAT");
                int action_type = kh_val(g_merge_config_token_2_idx, kitr);
                assert((action_type == MERGE_KEEP || action_type == MERGE_DROP || action_type == MERGE_MOVE_TO_FORMAT
                            || action_type == MERGE_COPY_TO_FORMAT)
                        && "Action to INFO/FORMAT must be KEEP|DROP|MOVE_TO_FORMAT|COPY_TO_FORMAT");
                if(token_type == MERGE_INFO)
                {
                    kh_put(strdict, g_merge_config.info_field2action, field, &ret);
                    kitr = kh_get(strdict, g_merge_config.info_field2action, field);
                    assert(kitr != kh_end(g_merge_config.info_field2action));
                    kh_val(g_merge_config.info_field2action, kitr) = action_type;
                }
                else
                {
                    kh_put(strdict, g_merge_config.format_field2action, field, &ret);
                    kitr = kh_get(strdict, g_merge_config.format_field2action, field);
                    assert(kitr != kh_end(g_merge_config.format_field2action));
                    kh_val(g_merge_config.format_field2action, kitr) = action_type;
                }
                free(ptr);
                break;
            default:
                fprintf(stderr,"Unexpected token %s at line %d .. skipping to next line\n",ptr,line_number);
                break;
        }
        ++line_number;
    }
    if(line != 0)
        free(line);
    fclose(fptr);
}

void destroy_merge_config()
{
    int i=0;
#define DESTROY_STRING_KEYS(h)                          \
    for(i=0;i<kh_end(h);++i)                            \
        if(kh_exist(h,i))                               \
            free((void*)kh_key(h,i));                
    if(g_merge_config_token_2_idx)
    {
        //strings are not dynamically allocated for this dict - const char* literals
        kh_destroy(strdict, g_merge_config_token_2_idx);
    }
    if(g_merge_config.info_field2action)
    {
        DESTROY_STRING_KEYS(g_merge_config.info_field2action);
        kh_destroy(strdict, g_merge_config.info_field2action);
    }
    if(g_merge_config.format_field2action)
    {
        DESTROY_STRING_KEYS(g_merge_config.format_field2action);
        kh_destroy(strdict, g_merge_config.format_field2action);
    }
    if(g_merge_config.file2id2action)
    {
        for(i=0;i<g_merge_config.num_files;++i)
            free(g_merge_config.file2id2action[i]);
        free(g_merge_config.file2id2action);
    }
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.\n");
    fprintf(stderr, "         Compatible records are combined into one according to the -m option.\n");
    fprintf(stderr, "Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        --force-samples                resolve duplicate sample names\n");
    fprintf(stderr, "        --print-header                 print only the merged header and exit\n");
    fprintf(stderr, "        --use-header <file>            use the provided header\n");
    fprintf(stderr, "    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or \"-\" to turn off the default [DP:sum,DP4:sum]\n");
    fprintf(stderr, "    -l, --file-list <file>             read file names from the file\n");
    fprintf(stderr, "    -m, --merge <string>               merge sites with differing alleles for <snps|indels|both|all|none|id>, see man page for details [both]\n");
    fprintf(stderr, "    -o, --output <file>                write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "\n");
    exit(1);
}

enum ArgsIdxEnum
{
  ARGS_IDX_MERGE_CONFIG_FILE=1000,
  ARGS_IDX_REFERENCE_FILE,
  ARGS_IDX_INPUT_GVCFS,
  ARGS_IDX_DO_GATK_MERGE,
  ARGS_IDX_MEASURE_ITERATOR_TIMING,
  ARGS_IDX_COMPUTE_PLMEDIAN,
  ARGS_IDX_TAG
};

int main_vcfmerge(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->collapse = COLLAPSE_BOTH;
    int regions_is_file = 0;
#ifdef DEBUG
    g_debug_fptr = fopen("debug.txt","w");
    g_vcf_debug_fptr = stderr;
    /*g_vcf_debug_fptr = g_debug_fptr;*/
    ks_resize(&g_debug_string, 4096);   //4KB
#endif

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"merge",1,0,'m'},
        {"file-list",1,0,'l'},
        {"apply-filters",1,0,'f'},
        {"use-header",1,0,1},
        {"print-header",0,0,2},
        {"force-samples",0,0,3},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"info-rules",1,0,'i'},
        {"merge-config",1,0,ARGS_IDX_MERGE_CONFIG_FILE},
        {"reference",1,0,ARGS_IDX_REFERENCE_FILE},
        {"gvcf",0,0,ARGS_IDX_INPUT_GVCFS},
        {"tag",1,0,ARGS_IDX_TAG},
        {"gatk",0,0,ARGS_IDX_DO_GATK_MERGE},
        {"iterator-timing",0,0,ARGS_IDX_MEASURE_ITERATOR_TIMING},
        {"PLmedian",1,0,ARGS_IDX_COMPUTE_PLMEDIAN},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hm:f:r:R:o:O:i:l:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'l': args->file_list = optarg; break;
            case 'i': args->info_rules = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'm':
                args->collapse = COLLAPSE_NONE;
                if ( !strcmp(optarg,"snps") ) args->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->collapse |= COLLAPSE_BOTH;
                else if ( !strcmp(optarg,"any") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"none") ) args->collapse = COLLAPSE_NONE;
                else if ( !strcmp(optarg,"id") ) { args->collapse = COLLAPSE_NONE; args->merge_by_id = 1; }
                else error("The -m type \"%s\" is not recognised.\n", optarg);
                break;
            case 'f': args->files->apply_filters = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case  1 : args->header_fname = optarg; break;
            case  2 : args->header_only = 1; break;
            case  3 : args->force_samples = 1; break;
            case ARGS_IDX_MERGE_CONFIG_FILE: parse_merge_config(optarg); break;
            case ARGS_IDX_REFERENCE_FILE: initialize_reference(args, optarg); break;
            case ARGS_IDX_INPUT_GVCFS: g_is_input_gvcf = 1; break;
            case ARGS_IDX_TAG:  args->m_tag = strdup(optarg); break;
            case ARGS_IDX_DO_GATK_MERGE: g_do_gatk_merge = 1; break;
            case ARGS_IDX_MEASURE_ITERATOR_TIMING: g_measure_iterator_timing_only = 1; break;
            case ARGS_IDX_COMPUTE_PLMEDIAN: args->m_plmedian_info.m_output_fptr = fopen(optarg, "w"); break;
            case 'h': 
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc==optind && !args->file_list ) usage();
    if ( argc-optind<2 && !args->file_list ) usage();

    assert((!g_is_input_gvcf || args->reference_filename) && "To merge gVCFs, you must provide a reference file using --reference=<filename.fa");

    args->files->require_index = 1;
    if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", args->regions_list);

    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open %s: %s\n", argv[optind],bcf_sr_strerror(args->files->errnum));
        optind++;
    }
    if ( args->file_list )
    {
        int nfiles, i;
        char **files = hts_readlines(args->file_list, &nfiles);
        if ( !files ) error("Failed to read from %s\n", args->file_list);
        for (i=0;i<nfiles; i++)
            if ( !bcf_sr_add_reader(args->files, files[i]) ) error("Failed to open %s: %s\n", files[i],bcf_sr_strerror(args->files->errnum));
        for (i=0; i<nfiles; i++) free(files[i]);
        free(files);
    }
#ifdef COLLECT_STATS
    initialize_stats(args);
#endif
    merge_vcf(args);
    //DESTRUCTION IS A FUNCTION!!!
    destroy_idmap(args);
    if(g_preprocess_vcfs)       //destroy modified headers
    {
        int i;
        for(i=0;i<args->files->nreaders;++i)
            bcf_hdr_destroy(args->files->readers[i].processed_header);
    }
    bcf_sr_destroy(args->files);
    destroy_merge_config();
    destroy_reference(args);
#ifdef COLLECT_STATS
    print_stats(args);
    destroy_stats(args);
#endif
    if(args->m_tag)
        free(args->m_tag);
    if(args->m_plmedian_info.m_output_fptr)
    {
        if(args->m_plmedian_info.m_buffer_len)
            free(args->m_plmedian_info.m_buffer);
        if(args->m_plmedian_info.m_median_result_len)
            free(args->m_plmedian_info.m_median_result);
        if(args->m_plmedian_info.m_reorg_buffer_len)
            free(args->m_plmedian_info.m_reorg_buffer);
        fclose(args->m_plmedian_info.m_output_fptr);
    }
    free(args);
#ifdef DEBUG
    fclose(g_debug_fptr);
    free(g_debug_string.s);
#endif
    return 0;
}

