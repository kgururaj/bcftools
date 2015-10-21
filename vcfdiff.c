#include <stdio.h>
#include <stdlib.h>
#include "htslib/vcf.h"
#include "htslib/khash.h"
#include <getopt.h>
#include <assert.h>
#include <math.h>

#ifdef DEBUG
#define ASSERT(X)  assert(X)
#else 
#define ASSERT(X)  ;
#endif

float g_FLOAT_TOLERANCE=1e-5f;

int g_printed_already = 0;
extern kstring_t g_debug_string;
int g_num_files = 0;
int g_skip_DP_info = 0;
int g_skip_AD_fmt = 0;
int g_skip_GT = 0;
typedef struct
{
    char* m_name;
    htsFile* m_fptr;
    bcf_hdr_t* m_hdr;
    bcf1_t* m_line;
}reader;
reader* g_files = 0;

typedef struct
{
    int* m_sample_map;
    int* m_id_map;
    int* m_contig_map;
    int* m_allele_map;
    int m_allele_map_size;
    int* m_gt_map;
    int m_gt_map_size;
    int m_GT_fmt_id;
    int m_AD_fmt_id;
    int m_DP_info_id;
}reader_map;

reader_map* g_readers_map;

int get_element_size(int bcf_ht_type);

//Check whether 2 files have different samples
void compare_sample_names(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    int j = 0;
    int k = 0;
    int first_is_first = (first_idx == 0);
    if(first_is_first && bcf_hdr_nsamples(second_file->m_hdr) != bcf_hdr_nsamples(first_file->m_hdr))
    {
        printf("ERROR: files %s and %s have different number of samples\n",first_file->m_name, second_file->m_name);
        exit(0);        //don't even bother doing anything else
    }
    for(j=0;j<second_file->m_hdr->n[BCF_DT_SAMPLE];++j)
    {
        const char* curr_sample_name = bcf_hdr_int2id(second_file->m_hdr, BCF_DT_SAMPLE, j);
        if(first_is_first)
            g_readers_map[second_idx].m_sample_map[j] = -1;
        int found_sample = 0;
        for(k=0;k<first_file->m_hdr->n[BCF_DT_SAMPLE];++k)
            if(strcmp(bcf_hdr_int2id(first_file->m_hdr, BCF_DT_SAMPLE, k), curr_sample_name) == 0)
            {
                found_sample = 1;
                if(first_is_first)
                    g_readers_map[second_idx].m_sample_map[j] = k;
                break;
            }
        if(!found_sample)
            printf("ERROR: file %s has sample %s which does not exist in file %s\n",second_file->m_name, curr_sample_name, first_file->m_name);
    }
}

//Check whether 2 files have different IDs
void compare_header_dictionary(reader* first_file, reader* second_file, int first_idx, int second_idx,
        int dict_idx)
{
    int j = 0;
    int k = 0;
    int first_is_first = (first_idx == 0);
    const char* dict_name = (dict_idx == BCF_DT_ID) ? "ID" : "contig";
    for(j=0;j<second_file->m_hdr->n[dict_idx];++j)
    {
        int* id_map = 0;
        if(first_is_first)
        {
            if(dict_idx == BCF_DT_ID)
                id_map = g_readers_map[second_idx].m_id_map;
            else
                id_map = g_readers_map[second_idx].m_contig_map;
            id_map[j] = -1;
        }
        const char* curr_key = bcf_hdr_int2id(second_file->m_hdr, dict_idx, j);
        if(curr_key == 0)
            continue;
        int found_id = 0;
        for(k=0;k<first_file->m_hdr->n[dict_idx];++k)
        {
            const char* first_file_key = bcf_hdr_int2id(first_file->m_hdr, dict_idx, k);
            if(first_file_key && strcmp(curr_key,first_file_key) == 0)
            {
                found_id = 1;
                if(first_is_first)
                    id_map[j] = k;
                break;
            }
        }
        if(!found_id)
            printf("NOTE: file %s has %s %s which does not exist in file %s\n",second_file->m_name, dict_name, curr_key, first_file->m_name);
    }
}


void print_lines(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    if(g_printed_already || first_idx != 0)
        return;
    g_printed_already = 1;
    printf("========================================================\n");
    vcf_format(first_file->m_hdr, first_file->m_line, &g_debug_string);
    printf("File %s : %s\n",first_file->m_name, g_debug_string.s);
    g_debug_string.l = 0;
    vcf_format(second_file->m_hdr, second_file->m_line, &g_debug_string);
    printf("File %s : %s\n",second_file->m_name, g_debug_string.s);
    g_debug_string.l = 0;
}

#define REALLOC_IF_NEEDED(ptr, curr_num_elements, new_num_elements, type)       \
    if((new_num_elements) > (curr_num_elements))                                \
    {                                                                           \
        (curr_num_elements) = (new_num_elements);                               \
        (ptr) = (type*)realloc((ptr), (new_num_elements)*sizeof(type));         \
    }

void compare_alleles(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    int j = 0;
    int k = 0;
    int first_is_first = (first_idx == 0);
    bcf1_t* first_line = first_file->m_line;
    bcf1_t* second_line = second_file->m_line;
    int diff_alleles = 0;
    if(first_is_first)
    {
        if(first_line->n_allele != second_line->n_allele)
        {
            print_lines(first_file, second_file, first_idx, second_idx);
            printf("Different alleles\n");
            return;
        }
        for(j=0;j<second_line->n_allele;++j)
        {
            char* second_allele = second_line->d.allele[j];
            int allele_found = 0;
            ASSERT(j < g_readers_map[second_idx].m_allele_map_size);
            for(k=0;k<first_line->n_allele;++k)
                if(strcmp(second_allele, first_line->d.allele[k]) == 0)
                {
                    allele_found = 1;
                    if(first_is_first)
                    {
                        g_readers_map[second_idx].m_allele_map[j] = k;
                    }
                    break;
                }
            if(!allele_found)
                diff_alleles = 1;
        }
        //GT map
        if(!diff_alleles)   //same alleles - setup genotype map
        {
            int orig_idx = 0;
            int map_idx = 0;
            for(j=0;j<second_line->n_allele;++j)
            {
                int map_j = g_readers_map[second_idx].m_allele_map[j];
                for(k=0;k<=j;++k)
                {
                    int map_k = g_readers_map[second_idx].m_allele_map[k];
                    orig_idx = bcf_alleles2gt(j, k);
                    map_idx = bcf_alleles2gt(map_j, map_k);
                    g_readers_map[second_idx].m_gt_map[orig_idx] = map_idx;
                }
            }
        }
        else
        {
            print_lines(first_file, second_file, first_idx, second_idx);
            printf("Different alleles\n");
        }
    }
}

int compare_unequal_char(char a, char b)
{
    return (a != b);
}

int compare_unequal_int(int32_t a, int32_t b)
{
    return (a != b);
}

int compare_unequal_float(float a, float b)
{
    if(bcf_float_is_missing(a) && bcf_float_is_missing(b))
        return 0;
    float abs_diff = fabsf(a-b);
    if(a != 0)
    {
        float rel_diff = fabsf(abs_diff/a);
        return (abs_diff > g_FLOAT_TOLERANCE && rel_diff > g_FLOAT_TOLERANCE);
    }
    else
        return (abs_diff > g_FLOAT_TOLERANCE);
}

char* g_diff_string = 0;
int compare_values(int coltype, reader* first_file, reader* second_file,
        int first_idx, int second_idx,
        int first_field_id, int second_field_id,
        void* first_values, int first_nvalues, void* second_values, int second_nvalues)
{
    if(first_nvalues != second_nvalues) //different lengths
    {
        g_diff_string = " different lengths ";
        return -1;
    }
    int first_type = bcf_hdr_id2type(first_file->m_hdr, coltype, first_field_id);
    int second_type = bcf_hdr_id2type(second_file->m_hdr, coltype, second_field_id);
    if(first_type != second_type)
    {
        g_diff_string = " different types ";
        return -1;
    }
    int first_length_type = bcf_hdr_id2length(first_file->m_hdr, coltype, first_field_id);
    int second_length_type = bcf_hdr_id2length(second_file->m_hdr, coltype, second_field_id);
    if(first_length_type != second_length_type)
    {
        g_diff_string = " different length types ";
        return -1;
    }
#define BRANCH(type_t, compare_unequal) \
    type_t* first_ptr = (type_t*)first_values; \
    type_t* second_ptr = (type_t*)second_values; \
    if(first_length_type == BCF_VL_A || first_length_type == BCF_VL_R) \
    { \
        int offset = (first_length_type == BCF_VL_A) ? 1 : 0; \
        for(i=0;i<second_nvalues;++i) \
        { \
            int mapped_idx = g_readers_map[second_idx].m_allele_map[i+offset] - offset; \
            if(mapped_idx < 0 || mapped_idx >= first_nvalues) \
            { \
                g_diff_string = " unmapped allele idx "; \
                return -1; \
            } \
            if(compare_unequal(first_ptr[mapped_idx], second_ptr[i])) \
            { \
                g_diff_string = " unequal values "; \
                return -1; \
            } \
        } \
    } \
    else \
        if(first_length_type == BCF_VL_G) \
        { \
            for(i=0;i<second_nvalues;++i) \
            { \
                int mapped_idx = g_readers_map[second_idx].m_gt_map[i]; \
                if(mapped_idx < 0 || mapped_idx >= first_nvalues) \
                { \
                    g_diff_string = " unmapped gt idx "; \
                    return -1; \
                } \
                if(compare_unequal(first_ptr[mapped_idx], second_ptr[i])) \
                { \
                    g_diff_string = " unequal values "; \
                    return -1; \
                } \
            } \
        } \
        else \
            for(i=0;i<second_nvalues;++i) \
            { \
                if(compare_unequal(first_ptr[i], second_ptr[i])) \
                { \
                    g_diff_string = " unequal values "; \
                    return -1; \
                } \
            }
    
    int i = 0;
    //Check for equality
    switch(first_type)
    {
        case BCF_HT_FLAG:
            return 1;   //both have same length, implies flag set
            break;
        case BCF_HT_INT:
            {
                BRANCH(int32_t, compare_unequal_int);
                return 1;
                break;
            }
        case BCF_HT_REAL:
            {
                BRANCH(float, compare_unequal_float);
                return 1;
                break;
            }
        case BCF_HT_STR:
            {
                BRANCH(char, compare_unequal_char);
                return 1;
                break;
            }
    }
    return -1;
}

void compare_info(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    int j = 0;
    int k = 0;
    int first_is_first = (first_idx == 0);
    bcf1_t* first_line = first_file->m_line;
    bcf1_t* second_line = second_file->m_line;
    if(first_is_first)
    {
        //Can't use n_info as some fields may be deleted
        int first_num_info_fields = 0;
        int second_num_info_fields = 0;
        for(j=0;j<first_line->n_info;++j)
        {
            if(first_line->d.info[j].vptr != 0)
                ++first_num_info_fields;
        }
        for(j=0;j<second_line->n_info;++j)
        {
            if(second_line->d.info[j].vptr != 0)
                ++second_num_info_fields;
        }
        if(first_num_info_fields != second_num_info_fields)
        {
            print_lines(first_file, second_file, first_idx, second_idx);
            printf("Different number of INFO fields\n");
            return;
        }
        
        char* first_values = 0;
        int first_allocated_nvalues = 0;
        int first_allocated_size = 0;
        int first_nvalues = 0;
        int first_element_size = 0;
        
        char* second_values = 0;
        int second_allocated_nvalues = 0;
        int second_allocated_size = 0;
        int second_nvalues = 0;
        int second_element_size = 0;
        for(j=0;j<second_line->n_info;++j)
        {
            bcf_info_t* second_info = &(second_line->d.info[j]);
            if(second_info->vptr == 0)
                continue;
	    if(g_skip_DP_info && second_info->key == g_readers_map[second_idx].m_DP_info_id)
		continue;
            second_element_size = get_element_size(bcf_hdr_id2type(second_file->m_hdr, BCF_HL_INFO, second_info->key));
            second_allocated_nvalues = second_allocated_size/second_element_size;
            second_nvalues = bcf_unpack_info_values(second_file->m_hdr, second_line, second_info,
                        (void**)&second_values, &second_allocated_nvalues,
                        bcf_hdr_id2type(second_file->m_hdr, BCF_HL_INFO, second_info->key));
            if(second_allocated_nvalues*second_element_size > second_allocated_size)
                second_allocated_size = second_allocated_nvalues*second_element_size;
            int mapped_info_id = g_readers_map[second_idx].m_id_map[second_info->key];
            int found_info = (mapped_info_id < 0) ? -1 : 0;
            if(found_info >= 0)
                for(k=0;k<first_line->n_info;++k)
                {
                    bcf_info_t* first_info = &(first_line->d.info[k]);
                    if(first_info->vptr == 0)
                        continue;
                    if(mapped_info_id == first_info->key)
                    {
                        first_element_size = get_element_size(bcf_hdr_id2type(first_file->m_hdr, BCF_HL_INFO, first_info->key));
                        first_allocated_nvalues = first_allocated_size/first_element_size;
                        first_nvalues = bcf_unpack_info_values(first_file->m_hdr, first_line, first_info,
                                (void**)&first_values, &first_allocated_nvalues,
                                bcf_hdr_id2type(first_file->m_hdr, BCF_HL_INFO, first_info->key));
                        if(first_allocated_nvalues*first_element_size > first_allocated_size)
                            first_allocated_size = first_allocated_nvalues*first_element_size;
                        int is_equal = compare_values(BCF_HL_INFO, first_file, second_file,
                                first_idx, second_idx, mapped_info_id, second_info->key,
                                (void*)first_values, first_nvalues, (void*)second_values, second_nvalues);
                        if(is_equal > 0)
                            found_info = 1;
                        break;
                    }
                }
            if(found_info <= 0)
            {
                print_lines(first_file, second_file, first_idx, second_idx);
                printf("Different INFO field %s\n",bcf_hdr_int2id(second_file->m_hdr, BCF_DT_ID, second_info->key));
            }
        }
        if(first_values)
            free(first_values);
        if(second_values)
            free(second_values);
    }
}

void compare_fmt(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    int j = 0;
    int k = 0;
    int l = 0;
    int first_is_first = (first_idx == 0);
    bcf1_t* first_line = first_file->m_line;
    bcf1_t* second_line = second_file->m_line;
    g_diff_string = 0;
    if(first_is_first)
    {
        //Can't use n_fmt as some fields may be deleted
        int first_num_fmt_fields = 0;
        int second_num_fmt_fields = 0;
        for(j=0;j<first_line->n_fmt;++j)
        {
            if(first_line->d.fmt[j].p != 0)
                ++first_num_fmt_fields;
        }
        for(j=0;j<second_line->n_fmt;++j)
        {
            if(second_line->d.fmt[j].p != 0)
                ++second_num_fmt_fields;
        }
        if(first_num_fmt_fields != second_num_fmt_fields)
        {
            print_lines(first_file, second_file, first_idx, second_idx);
            printf("Different number of FORMAT fields\n");
            return;
        }
        
        char* first_values = 0;
        int first_allocated_nvalues = 0;
        int first_allocated_size = 0;
        int first_nvalues = 0;
        int first_element_size = 0;
        int first_per_sample_nvalues = 0;
        int first_offset = 0;
        
        char* second_values = 0;
        int second_allocated_nvalues = 0;
        int second_allocated_size = 0;
        int second_nvalues = 0;
        int second_element_size = 0;
        int second_per_sample_nvalues = 0;
        int second_offset = 0;
        ASSERT(bcf_hdr_nsamples(first_file->m_hdr) == bcf_hdr_nsamples(second_file->m_hdr));
        int nsamples = bcf_hdr_nsamples(first_file->m_hdr);
        for(j=0;j<second_line->n_fmt;++j)
        {
            bcf_fmt_t* second_fmt = &(second_line->d.fmt[j]);
            if(second_fmt->p == 0 || (g_skip_GT && second_fmt->id == g_readers_map[second_idx].m_GT_fmt_id)
		|| (g_skip_AD_fmt && second_fmt->id == g_readers_map[second_idx].m_AD_fmt_id) )
                continue;
            second_element_size = get_element_size(bcf_hdr_id2type(second_file->m_hdr, BCF_HL_FMT, second_fmt->id));
            second_allocated_nvalues = second_allocated_size/second_element_size;
            second_nvalues = bcf_unpack_format_values(second_file->m_hdr, second_line, second_fmt,
                        (void**)&second_values, &second_allocated_nvalues,
                        bcf_hdr_id2type(second_file->m_hdr, BCF_HL_FMT, second_fmt->id));
            if(second_allocated_nvalues*second_element_size > second_allocated_size)
                second_allocated_size = second_allocated_nvalues*second_element_size;
            second_per_sample_nvalues = second_nvalues/nsamples;
            int mapped_fmt_id = g_readers_map[second_idx].m_id_map[second_fmt->id];
            int found_fmt = (mapped_fmt_id < 0) ? -1 : 0;
            if(found_fmt >= 0)
                for(k=0;k<first_line->n_fmt;++k)
                {
                    bcf_fmt_t* first_fmt = &(first_line->d.fmt[k]);
                    if(first_fmt->p == 0)
                        continue;
                    if(mapped_fmt_id == first_fmt->id)
                    {
                        first_element_size = get_element_size(bcf_hdr_id2type(first_file->m_hdr, BCF_HL_FMT, first_fmt->id));
                        first_allocated_nvalues = first_allocated_size/first_element_size;
                        first_nvalues = bcf_unpack_format_values(first_file->m_hdr, first_line, first_fmt,
                                (void**)&first_values, &first_allocated_nvalues,
                                bcf_hdr_id2type(first_file->m_hdr, BCF_HL_FMT, first_fmt->id));
                        if(first_allocated_nvalues*first_element_size > first_allocated_size)
                            first_allocated_size = first_allocated_nvalues*first_element_size;
                        first_per_sample_nvalues = first_nvalues/nsamples;
                        found_fmt = 1;
                        //Compare one sample at a time
                        for(l=0;l<nsamples;++l)
                        {
                            int mapped_sample_idx = g_readers_map[second_idx].m_sample_map[l];
                            first_offset = mapped_sample_idx*first_per_sample_nvalues*first_element_size;
                            second_offset = l*second_per_sample_nvalues*second_element_size;
                            ASSERT(mapped_sample_idx >= 0 && mapped_sample_idx < nsamples);
                            int is_equal = compare_values(BCF_HL_FMT, first_file, second_file,
                                    first_idx, second_idx, mapped_fmt_id, second_fmt->id,
                                    (void*)(first_values+first_offset), first_per_sample_nvalues,
                                    (void*)(second_values+second_offset), second_per_sample_nvalues);
                            if(is_equal < 0)
                            {
                                found_fmt = 0;
                                break;
                            }
                        }
                        break;
                    }
                }
            else
                g_diff_string = " invalid map to golden ";
            if(found_fmt <= 0)
            {
                print_lines(first_file, second_file, first_idx, second_idx);
                printf("Different FORMAT field %s diff string : %s\n",bcf_hdr_int2id(second_file->m_hdr, BCF_DT_ID, second_fmt->id), g_diff_string);
            }
        }
        if(first_values)
            free(first_values);
        if(second_values)
            free(second_values);
    }
}

void compare_filters(reader* first_file, reader* second_file, int first_idx, int second_idx)
{
    int j = 0;
    int k = 0;
    int first_is_first = (first_idx == 0);
    bcf1_t* first_line = first_file->m_line;
    bcf1_t* second_line = second_file->m_line;
    if(first_line->d.n_flt != second_line->d.n_flt)
    {
        print_lines(first_file, second_file, first_idx, second_idx);
        printf("Different FILTER fields\n");
        return;
    }
    if(first_is_first)
    {
        int diff_filter = 0;
        for(j=0;j<second_line->d.n_flt;++j)
        {
            ASSERT(second_line->d.flt[j] < second_file->m_hdr->n[BCF_DT_ID]);
            int mapped_flt_id = g_readers_map[second_idx].m_id_map[second_line->d.flt[j]];
            int found_flt = (mapped_flt_id < 0) ? -1 : 0;
            if(found_flt >= 0)
                for(k=0;k<first_line->d.n_flt;++k)
                    if(mapped_flt_id == first_line->d.flt[k])
                    {
                        found_flt = 1;
                        break;
                    }
            if(found_flt <= 0)
                diff_filter = 1;
        }
        if(diff_filter)
        {
            print_lines(first_file, second_file, first_idx, second_idx);
            printf("Different FILTER fields\n");
        }
    }
}

enum ArgsIdxEnum
{
  ARGS_IDX_SKIP_DP_INFO=1000,
  ARGS_IDX_SKIP_AD_FMT,
  ARGS_IDX_SKIP_GT
};

int main_vcfdiff(int argc, char *argv[])
{
    int c ;
    static struct option loptions[] =
    {
        {"skip-DP-info",0,0,ARGS_IDX_SKIP_DP_INFO},
        {"skip-AD-fmt",0,0,ARGS_IDX_SKIP_AD_FMT},
        {"skip-GT",0,0,ARGS_IDX_SKIP_GT},
        {"tolerance",0,0,'t'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "t:",loptions,NULL)) >= 0) {
	switch(c)
	{
            case 't':
                g_FLOAT_TOLERANCE = strtof(optarg, 0);
                break;
	    case ARGS_IDX_SKIP_DP_INFO:
		g_skip_DP_info = 1;
		break;
	    case ARGS_IDX_SKIP_AD_FMT:
		g_skip_AD_fmt = 1;
		break;
            case ARGS_IDX_SKIP_GT:
                g_skip_GT = 1;
                break;
	    default:
		fprintf(stderr,"Unknown option\n");
		exit(-1);
	}
    }
    assert(argc - optind >= 2 && "Requires at least 2 arguments <vcf1.gz> <vcf2.gz>");
    argv = argv + optind;
    g_num_files = argc - optind;
    g_files = (reader*)malloc(g_num_files*sizeof(reader));
    g_readers_map = (reader_map*)calloc(g_num_files, sizeof(reader_map));
    //debug string
    ks_resize(&g_debug_string, 4096);   //4KB
    //open file and read record
    int i = 0;
    reader* first_file = &(g_files[0]);
    for(i=0;i<g_num_files;++i)
    {
        reader* curr_file = &(g_files[i]);
        curr_file->m_fptr = hts_open(argv[i], "r");
        curr_file->m_name = argv[i];
        assert(curr_file->m_fptr);
        curr_file->m_hdr = bcf_hdr_read(curr_file->m_fptr);
        curr_file->m_line = bcf_init();
        assert(curr_file->m_line);
        //sample map
        g_readers_map[i].m_sample_map = (int*)malloc(bcf_hdr_nsamples(curr_file->m_hdr)*sizeof(int));
        compare_sample_names(first_file, curr_file, 0, i);
        compare_sample_names(curr_file, first_file, i, 0);
        
        //idmap
        g_readers_map[i].m_id_map = (int*)malloc(curr_file->m_hdr->n[BCF_DT_ID]*sizeof(int));
        compare_header_dictionary(first_file, curr_file, 0, i, BCF_DT_ID);
        compare_header_dictionary(curr_file, first_file, i, 0, BCF_DT_ID);

        //contig map        
        g_readers_map[i].m_contig_map = (int*)malloc(curr_file->m_hdr->n[BCF_DT_CTG]*sizeof(int));
        compare_header_dictionary(first_file, curr_file, 0, i, BCF_DT_CTG);
        compare_header_dictionary(curr_file, first_file, i, 0, BCF_DT_CTG);

        //allele map space allocation
        g_readers_map[i].m_allele_map_size = 10;
        g_readers_map[i].m_allele_map = (int*)malloc(g_readers_map[i].m_allele_map_size*sizeof(int));
        //genotype map
        g_readers_map[i].m_gt_map_size = 
            (g_readers_map[i].m_allele_map_size*(g_readers_map[i].m_allele_map_size+1))/2;
        g_readers_map[i].m_gt_map = (int*)malloc(g_readers_map[i].m_gt_map_size*sizeof(int));
        //GT id
        g_readers_map[i].m_GT_fmt_id = bcf_hdr_id2int(curr_file->m_hdr, BCF_DT_ID, "GT");
        ASSERT(g_readers_map[i].m_GT_fmt_id);
        //AD id
        g_readers_map[i].m_AD_fmt_id = bcf_hdr_id2int(curr_file->m_hdr, BCF_DT_ID, "AD");
        ASSERT(g_readers_map[i].m_AD_fmt_id);
	//DP id
        g_readers_map[i].m_DP_info_id = bcf_hdr_id2int(curr_file->m_hdr, BCF_DT_ID, "DP");
        ASSERT(g_readers_map[i].m_DP_info_id);
        //set length type of AD to R
        int AD_id = bcf_hdr_id2int(curr_file->m_hdr, BCF_DT_ID, "AD");
        if(AD_id >= 0)
            bcf_hdr_set_length_descriptor(curr_file->m_hdr, BCF_HL_FMT, AD_id, BCF_VL_R);
    }
    int first_status = 0;
    while(first_status >= 0)
    {
        bcf_clear(first_file->m_line);
        first_status = bcf_read(first_file->m_fptr, first_file->m_hdr, first_file->m_line);
        bcf1_t* first_line = first_file->m_line;
        if(first_status >= 0)
            bcf_unpack(first_file->m_line, BCF_UN_ALL);
        for(i=1;i<g_num_files;++i)
        {
            reader* curr_file = &(g_files[i]);
            reader_map* curr_map = &(g_readers_map[i]);
            bcf_clear(curr_file->m_line);
            int curr_status = bcf_read(curr_file->m_fptr, curr_file->m_hdr, curr_file->m_line);
            if(curr_status >= 0)
                bcf_unpack(curr_file->m_line, BCF_UN_ALL);
            bcf1_t* curr_line = curr_file->m_line;
            if((curr_status < 0 && first_status >= 0) || (first_status < 0 && curr_status >= 0))
            {
                printf("Different number of lines in file %s and %s\n", first_file->m_name, curr_file->m_name);
                continue;
            }
            if(curr_status < 0)
                continue;
            //reset maps
            g_printed_already = 0;
            int num_alleles = curr_file->m_line->n_allele;
            REALLOC_IF_NEEDED(curr_map->m_allele_map, curr_map->m_allele_map_size, num_alleles, int);
            memset(curr_map->m_allele_map, -1, curr_map->m_allele_map_size*sizeof(int));
            int num_gts = (num_alleles*(num_alleles+1))/2;
            REALLOC_IF_NEEDED(curr_map->m_gt_map, curr_map->m_gt_map_size, num_gts, int);
            memset(curr_map->m_gt_map, -1, curr_map->m_gt_map_size*sizeof(int));
            //chrom
            if(curr_map->m_contig_map[curr_line->rid] != first_line->rid)
            {
                print_lines(first_file, curr_file, 0, i);
                printf("Different CHROM values\n");
            }
            //pos
            if(curr_line->pos != first_line->pos)
            {
                print_lines(first_file, curr_file, 0, i);
                printf("Different POS values\n");
            }
            //ignoring ID column
            //qual
            if(!(bcf_float_is_missing(curr_line->qual) && bcf_float_is_missing(first_line->qual))
                    && curr_line->qual != first_line->qual)
            {
                print_lines(first_file, curr_file, 0, i);
                printf("Different QUAL values %f %f\n",curr_line->qual, first_line->qual);
            }
            compare_alleles(first_file, curr_file, 0, i);
            compare_filters(first_file, curr_file, 0, i);
            compare_info(first_file, curr_file, 0, i);
            compare_fmt(first_file, curr_file, 0, i);
            if(g_printed_already)
                printf("========================================================\n");
        }
    }
   
    //Free! 
    for(i=0;i<g_num_files;++i)
    {
        free(g_readers_map[i].m_id_map);
        free(g_readers_map[i].m_contig_map);
        free(g_readers_map[i].m_sample_map);
        free(g_readers_map[i].m_allele_map);
        free(g_readers_map[i].m_gt_map);

        bcf_destroy(g_files[i].m_line);
        bcf_hdr_destroy(g_files[i].m_hdr);
        hts_close(g_files[i].m_fptr);
    }
    free(g_debug_string.s);
    free(g_readers_map);
    free(g_files);
    return 0;
}
