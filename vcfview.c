/*  vcfview.c -- VCF/BCF conversion, view, subset and filter VCF/BCF files.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

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
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "htslib/khash_str2int.h"
#include <inttypes.h>
#include <time.h>
#include <float.h>

#ifdef DEBUG
#define ASSERT(x)  assert(x)
#else
#define ASSERT(x) ;
#endif

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define ALLELE_NONREF 1
#define ALLELE_MINOR 2
#define ALLELE_ALT1 3
#define ALLELE_MAJOR 4
#define ALLELE_NONMAJOR 5

#define GT_NEED_HOM 1
#define GT_NEED_HET 2
#define GT_NO_HOM   3
#define GT_NO_HET   4
#define GT_NEED_MISSING 5
#define GT_NO_MISSING 6

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_list, *regions_list;
    int argc, clevel, output_type, print_header, update_info, header_only, n_samples, *imap, calc_ac;
    int trim_alts, sites_only, known, novel, min_alleles, max_alleles, private_vars, uncalled, phased;
    int min_ac, min_ac_type, max_ac, max_ac_type, min_af_type, max_af_type, gt_type;
    int *ac, mac;
    float min_af, max_af;
    char *fn_ref, *fn_out, **samples;
    int sample_is_file, force_samples;
    char *include_types, *exclude_types;
    int include, exclude;
    htsFile *out;
    sqlite_mappings_struct m_mapping_info;
    plmedian_struct  m_plmedian_info;   //0 in calloc
    csv_output_struct m_csv_output_info;
    int m_first_iteration;
    char* m_query_positions_file;
    char* m_apply_filters;
}
args_t;

void* allocate_sqlite3_mapping(const char* sqlite_file)
{
    sqlite_mappings_struct* mapping_info = (sqlite_mappings_struct*)calloc(1, sizeof(sqlite_mappings_struct));
    assert(mapping_info);
    strcpy(mapping_info->sqlite_file, sqlite_file);
    open_sqlite3_db(sqlite_file, &(mapping_info->db));
    return (void*)mapping_info;
}

void open_sqlite3_db(const char* sqlite_file, sqlite3** db)
{
    if(sqlite_file[0] == '\0')
    {
	fprintf(stderr, "To print out a TileDB csv file, you must pass the path to the samples and field names sqlite table using the --sqlite argument\n");
	exit(-1);
    }
    if(sqlite3_open(sqlite_file, db) != 0)
    {
	fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(*db)); 
	sqlite3_close(*db);
	exit(-1);
    }
}

static void init_data(args_t *args)
{
    int i;
    args->hdr = args->files->readers[0].header;

    if (args->calc_ac && args->update_info)
    {
        bcf_hdr_append(args->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
        bcf_hdr_append(args->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    }
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_view");

    if(args->m_first_iteration == 0)
	return;

    // setup sample data
    if (args->sample_names)
    {
        void *hdr_samples = khash_str2int_init();
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
            khash_str2int_inc(hdr_samples, bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i));

        void *exclude = (args->sample_names[0]=='^') ? khash_str2int_init() : NULL;
        int nsmpl;
        char **smpl = NULL;
        args->samples = NULL; args->n_samples = 0;
        smpl = hts_readlist(exclude ? &args->sample_names[1] : args->sample_names, args->sample_is_file, &nsmpl);
        if ( !smpl )
        {
            error("Could not read the list: \"%s\"\n", exclude ? &args->sample_names[1] : args->sample_names);
        }

        if ( exclude )
        {
            for (i=0; i<nsmpl; i++) {
                if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                    if (args->force_samples) {
                        fprintf(stderr, "Warn: exclude called for sample that does not exist in header: \"%s\"... skipping\n", smpl[i]);
                    } else {
                        error("Error: exclude called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                }
                khash_str2int_inc(exclude, smpl[i]);
            }

            for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
            {
                if ( exclude && khash_str2int_has_key(exclude,bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i))  ) continue;
                args->samples = (char**) realloc(args->samples, (args->n_samples+1)*sizeof(const char*));
                args->samples[args->n_samples++] = strdup(bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i));
            }
            khash_str2int_destroy(exclude);
        }
        else
        {
            for (i=0; i<nsmpl; i++) {
                if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                    if (args->force_samples) {
                        fprintf(stderr, "Warn: subset called for sample that does not exist in header: \"%s\"... skipping\n", smpl[i]);
                        continue;
                    } else {
                        error("Error: subset called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                }
                args->samples = (char**) realloc(args->samples, (args->n_samples+1)*sizeof(const char*));
                args->samples[args->n_samples++] = strdup(smpl[i]);
            }
        }
        for (i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
        khash_str2int_destroy(hdr_samples);
        if (args->n_samples == 0) {
            fprintf(stderr, "Warn: subsetting has removed all samples\n");
            args->sites_only = 1;
        }
    }

    if (args->n_samples)
        args->imap = (int*)malloc(args->n_samples * sizeof(int));

    // determine variant types to include/exclude
    if (args->include_types || args->exclude_types) {
        if (args->include_types && args->exclude_types) {
            fprintf(stderr, "Error: only supply one of --include-types, --exclude-types options\n");
            exit(1);
        }
        char **type_list = 0;
        int m = 0, n = 0;
        const char *q, *p;
        for (q = p = args->include_types ? args->include_types : args->exclude_types;; ++p) {
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    type_list = (char**)realloc(type_list, m * sizeof(char*));
                }
                type_list[n] = (char*)calloc(p - q + 1, 1);
                strncpy(type_list[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
        }
        type_list = (char**)realloc(type_list, n * sizeof(char*));

        if (args->include_types) {
            args->include = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->include |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->include |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->include |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->include |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    fprintf(stderr, "Accepted types are snps, indels, mnps, other\n");
                    exit(1);
                }
            }
        }
        if (args->exclude_types) {
            args->exclude = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->exclude |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->exclude |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->exclude |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->exclude |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    fprintf(stderr, "Accepted types are snps, indels, mnps, other\n");
                    exit(1);
                }
            }
        }
        for (i = 0; i < n; ++i)
            free(type_list[i]);
        free(type_list);
    }

    // setup output
    char modew[8];
    strcpy(modew, "w");
    if(args->output_type != FT_TILEDB_CSV)
    {
        if (args->clevel >= 0 && args->clevel <= 9) sprintf(modew + 1, "%d", args->clevel);
        if (args->output_type==FT_BCF) strcat(modew, "bu");         // uncompressed BCF
        else if (args->output_type & FT_BCF) strcat(modew, "b");    // compressed BCF
        else if (args->output_type & FT_GZ) strcat(modew,"z");      // compressed VCF
        args->out = hts_open(args->fn_out ? args->fn_out : "-", modew);
        if ( !args->out ) error("%s: %s\n", args->fn_out,strerror(errno));
    }
    else
    {
	if(args->m_mapping_info.db == 0)
	    open_sqlite3_db(args->m_mapping_info.sqlite_file, &(args->m_mapping_info.db));
    }

    // headers: hdr=full header, hsub=subset header, hnull=sites only header
    if (args->sites_only){
        args->hnull = bcf_hdr_subset(args->hdr, 0, 0, 0);
        bcf_hdr_remove(args->hnull, BCF_HL_FMT, NULL);
    }
    if (args->n_samples > 0)
    {
        args->hsub = bcf_hdr_subset(args->hdr, args->n_samples, args->samples, args->imap);
        if ( !args->hsub ) error("Error occurred while subsetting samples\n");
        if ( args->n_samples != bcf_hdr_nsamples(args->hsub) )
        {
            int i;
            for (i=0; i<args->n_samples; i++)
                if ( args->imap[i]<0 ) error("Error: No such sample: \"%s\"\n", args->samples[i]);
        }
    }

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}

void initialize_csv_output_info(csv_output_struct* ptr, FILE* output_fptr, unsigned buffer_size, uint8_t reset)
{
    if(reset)
        memset(ptr, 0, sizeof(csv_output_struct));
    ptr->m_csv_out_fptr = output_fptr;
    ptr->m_csv_out_buffer.m_buffer = (char*)malloc(buffer_size);
    ptr->m_csv_out_buffer.m_size = buffer_size;
    ptr->m_csv_out_buffer.m_offset = 0;
    ptr->m_htslib_buffer_size = 16384;
    ptr->m_htslib_buffer = (char*)malloc(16384);
}

void free_sqlite3_data(void* info_ptr)
{
    sqlite_mappings_struct* mapping_info = (sqlite_mappings_struct*)info_ptr;
    if(mapping_info->input_field_idx_2_global_idx)
	free(mapping_info->input_field_idx_2_global_idx);
    if(mapping_info->input_sample_idx_2_global_idx)
	free(mapping_info->input_sample_idx_2_global_idx);
    if(mapping_info->input_contig_idx_2_global_idx)
	free(mapping_info->input_contig_idx_2_global_idx);
    if(mapping_info->input_contig_idx_2_offset)
	free(mapping_info->input_contig_idx_2_offset);
    if(mapping_info->db)
	sqlite3_close(mapping_info->db);
    int64_t i = 0;
    if(mapping_info->m_field_names)
    {
        if(mapping_info->m_field_name_strings_allocated)
            for(i=0;i<mapping_info->m_num_fields;++i)
                if(mapping_info->m_field_names[i])
                    free(mapping_info->m_field_names[i]);
        free(mapping_info->m_field_names);
    }
    if(mapping_info->m_contig_names)
    {
        if(mapping_info->m_contig_name_strings_allocated)
            for(i=0;i<mapping_info->m_num_contigs;++i)
                if(mapping_info->m_contig_names[i])
                    free(mapping_info->m_contig_names[i]);
        free(mapping_info->m_contig_names);
    }
    if(mapping_info->m_contig_lengths)
        free(mapping_info->m_contig_lengths);
    if(mapping_info->m_sample_names)
    {
        if(mapping_info->m_sample_name_strings_allocated)
            for(i=0;i<mapping_info->m_num_samples;++i)
                if(mapping_info->m_sample_names[i])
                    free(mapping_info->m_sample_names[i]);
        free(mapping_info->m_sample_names);
    }
    if(mapping_info->m_tiledb_override_sample_name)
        free(mapping_info->m_tiledb_override_sample_name);
    memset(mapping_info, 0, sizeof(sqlite_mappings_struct));
}

void free_csv_output_info(csv_output_struct* ptr)
{
    if(ptr->m_csv_out_fptr)
    {
        fflush(ptr->m_csv_out_fptr);
        fclose(ptr->m_csv_out_fptr);
    }
    if(ptr->m_csv_out_buffer.m_buffer)
      free(ptr->m_csv_out_buffer.m_buffer);
    if(ptr->m_htslib_buffer)
        free(ptr->m_htslib_buffer);
    memset(ptr, 0, sizeof(csv_output_struct));
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->imap ) {
        for (i = 0; i < args->n_samples; ++i)
            free(args->samples[i]);
        free(args->samples);
        free(args->imap);
    }
    if (args->hnull) bcf_hdr_destroy(args->hnull);
    if (args->hsub) bcf_hdr_destroy(args->hsub);
    if ( args->filter )
        filter_destroy(args->filter);
    free(args->ac);
    free_sqlite3_data(&(args->m_mapping_info));
    free_csv_output_info(&(args->m_csv_output_info));
}

// true if all samples are phased.
// haploid genotypes are considered phased
// ./. => not phased, .|. => phased
int bcf_all_phased(const bcf_hdr_t *header, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_ptr = bcf_get_fmt(header, line, "GT");
    int all_phased = 1;
    if ( fmt_ptr )
    {
        int i, isample;
        for (isample=0; isample<line->n_sample; isample++)
        {
            int sample_phased = 0;
            #define BRANCH_INT(type_t,vector_end) { \
                type_t *p = (type_t*) (fmt_ptr->p + isample*fmt_ptr->size); \
                for (i=0; i<fmt_ptr->n; i++) \
                { \
                    if (fmt_ptr->n == 1 || (p[i] == vector_end && i == 1)) { sample_phased = 1; break; } /* haploid phased by definition */ \
                    if ( p[i] == vector_end ) { break; }; /* smaller ploidy */ \
                    if ( bcf_gt_is_missing(p[i]) ) continue; /* missing allele */ \
                    if ((p[i])&1) { \
                        sample_phased = 1; \
                        break; \
                    } \
                } \
            }
            switch (fmt_ptr->type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
                default: fprintf(stderr, "[E::%s] todo: fmt_type %d\n", __func__, fmt_ptr->type); exit(1); break;
            }
            #undef BRANCH_INT
            if (!sample_phased) {
                all_phased = 0;
                break;
            }
        }
    }
    return all_phased;
}

int subset_vcf(args_t *args, bcf1_t *line)
{
    if ( args->min_alleles && line->n_allele < args->min_alleles ) return 0; // min alleles
    if ( args->max_alleles && line->n_allele > args->max_alleles ) return 0; // max alleles
    if (args->novel || args->known)
    {
        if ( args->novel && (line->d.id[0]!='.' || line->d.id[1]!=0) ) return 0; // skip sites which are known, ID != '.'
        if ( args->known && line->d.id[0]=='.' && line->d.id[1]==0 ) return 0;  // skip sites which are novel, ID == '.'
    }

    if (args->include || args->exclude)
    {
        int line_type = bcf_get_variant_types(line);
        if ( args->include && !(line_type&args->include) ) return 0; // include only given variant types
        if ( args->exclude &&   line_type&args->exclude  ) return 0; // exclude given variant types
    }

    if ( args->filter )
    {
        int ret = filter_test(args->filter, line, NULL);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) return 0; }
        else if ( ret ) return 0;
    }

    hts_expand(int, line->n_allele, args->mac, args->ac);
    int i, an = 0, non_ref_ac = 0;
    if (args->calc_ac) {
        bcf_calc_ac(args->hdr, line, args->ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
        for (i=1; i<line->n_allele; i++)
            non_ref_ac += args->ac[i];
        for (i=0; i<line->n_allele; i++)
            an += args->ac[i];
    }

    if (args->n_samples)
    {
        int non_ref_ac_sub = 0, *ac_sub = (int*) calloc(line->n_allele,sizeof(int));
        bcf_subset(args->hdr, line, args->n_samples, args->imap);
        if (args->calc_ac) {
            bcf_calc_ac(args->hsub, line, ac_sub, BCF_UN_FMT); // recalculate AC and AN
            an = 0;
            for (i=0; i<line->n_allele; i++) {
                args->ac[i] = ac_sub[i];
                an += ac_sub[i];
            }
            for (i=1; i<line->n_allele; i++)
                non_ref_ac_sub += ac_sub[i];
            if (args->private_vars) {
                if (args->private_vars == FLT_INCLUDE && !(non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub)) { free(ac_sub); return 0; } // select private sites
                if (args->private_vars == FLT_EXCLUDE && non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub) { free(ac_sub); return 0; } // exclude private sites
            }
            non_ref_ac = non_ref_ac_sub;
        }
        free(ac_sub);
    }

    bcf_fmt_t *gt_fmt;
    if ( args->gt_type && (gt_fmt=bcf_get_fmt(args->hdr,line,"GT")) )
    {
        int nhet = 0, nhom = 0, nmiss = 0;
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
        {
            int type = bcf_gt_type(gt_fmt,i,NULL,NULL);
            if ( type==GT_HET_RA || type==GT_HET_AA )
            {
                if ( args->gt_type==GT_NO_HET ) return 0;
                nhet = 1;
            }
            else if ( type==GT_UNKN )
            {
                if ( args->gt_type==GT_NO_MISSING ) return 0;
                nmiss = 1;
            }
            else
            {
                if ( args->gt_type==GT_NO_HOM ) return 0;
                nhom = 1;
            }
        }
        if ( args->gt_type==GT_NEED_HOM && !nhom ) return 0;
        else if ( args->gt_type==GT_NEED_HET && !nhet ) return 0;
        else if ( args->gt_type==GT_NEED_MISSING && !nmiss ) return 0;
    }

    int minor_ac = 0;
    int major_ac = 0;
    if ( args->calc_ac )
    {
        minor_ac = args->ac[0];
        major_ac = args->ac[0];
        for (i=1; i<line->n_allele; i++){
            if (args->ac[i] < minor_ac) { minor_ac = args->ac[i]; }
            if (args->ac[i] > major_ac) { major_ac = args->ac[i]; }
        }
    }

    if (args->min_ac)
    {
        if (args->min_ac_type == ALLELE_NONREF && args->min_ac>non_ref_ac) return 0; // min AC
        else if (args->min_ac_type == ALLELE_MINOR && args->min_ac>minor_ac) return 0; // min minor AC
        else if (args->min_ac_type == ALLELE_ALT1 && args->min_ac>args->ac[1]) return 0; // min 1st alternate AC
        else if (args->min_ac_type == ALLELE_MAJOR && args->min_ac > major_ac) return 0; // min major AC
        else if (args->min_ac_type == ALLELE_NONMAJOR && args->min_ac > an-major_ac) return 0; // min non-major AC
    }
    if (args->max_ac)
    {
        if (args->max_ac_type == ALLELE_NONREF && args->max_ac<non_ref_ac) return 0; // max AC
        else if (args->max_ac_type == ALLELE_MINOR && args->max_ac<minor_ac) return 0; // max minor AC
        else if (args->max_ac_type == ALLELE_ALT1 && args->max_ac<args->ac[1]) return 0; // max 1st alternate AC
        else if (args->max_ac_type == ALLELE_MAJOR && args->max_ac < major_ac) return 0; // max major AC
        else if (args->max_ac_type == ALLELE_NONMAJOR && args->max_ac < an-major_ac) return 0; // max non-major AC
    }
    if (args->min_af)
    {
        if (an == 0) return 0; // freq not defined, skip site
        if (args->min_af_type == ALLELE_NONREF && args->min_af>non_ref_ac/(double)an) return 0; // min AF
        else if (args->min_af_type == ALLELE_MINOR && args->min_af>minor_ac/(double)an) return 0; // min minor AF
        else if (args->min_af_type == ALLELE_ALT1 && args->min_af>args->ac[1]/(double)an) return 0; // min 1st alternate AF
        else if (args->min_af_type == ALLELE_MAJOR && args->min_af > major_ac/(double)an) return 0; // min major AF
        else if (args->min_af_type == ALLELE_NONMAJOR && args->min_af > (an-major_ac)/(double)an) return 0; // min non-major AF
    }
    if (args->max_af)
    {
        if (an == 0) return 0; // freq not defined, skip site
        if (args->max_af_type == ALLELE_NONREF && args->max_af<non_ref_ac/(double)an) return 0; // max AF
        else if (args->max_af_type == ALLELE_MINOR && args->max_af<minor_ac/(double)an) return 0; // max minor AF
        else if (args->max_af_type == ALLELE_ALT1 && args->max_af<args->ac[1]/(double)an) return 0; // max 1st alternate AF
        else if (args->max_af_type == ALLELE_MAJOR && args->max_af < major_ac/(double)an) return 0; // max major AF
        else if (args->max_af_type == ALLELE_NONMAJOR && args->max_af < (an-major_ac)/(double)an) return 0; // max non-major AF
    }
    if (args->uncalled) {
        if (args->uncalled == FLT_INCLUDE && an > 0) return 0; // select uncalled
        if (args->uncalled == FLT_EXCLUDE && an == 0) return 0; // skip if uncalled
    }
    if (args->calc_ac && args->update_info) {
        bcf_update_info_int32(args->hdr, line, "AC", &args->ac[1], line->n_allele-1);
        bcf_update_info_int32(args->hdr, line, "AN", &an, 1);
    }
    if (args->trim_alts)
    {
        int ret = bcf_trim_alleles(args->hsub ? args->hsub : args->hdr, line);
        if ( ret==-1 ) error("Error: some GT index is out of bounds at %s:%d\n", bcf_seqname(args->hsub ? args->hsub : args->hdr, line), line->pos+1);
    }
    if (args->phased) {
        int phased = bcf_all_phased(args->hdr, line);
        if (args->phased == FLT_INCLUDE && !phased) { return 0; } // skip unphased
        if (args->phased == FLT_EXCLUDE && phased) { return 0; } // skip phased
    }
    if (args->sites_only) bcf_subset(args->hsub ? args->hsub : args->hdr, line, 0, 0);
    return 1;
}

void set_allele_type (int *atype, char *atype_string)
{
    *atype = ALLELE_NONREF;
    if (strcmp(atype_string, "minor") == 0) {
        *atype = ALLELE_MINOR;
    }
    else if (strcmp(atype_string, "alt1") == 0) {
        *atype = ALLELE_ALT1;
    }
    else if (strcmp(atype_string, "nref") == 0) {
        *atype = ALLELE_NONREF;
    }
    else if (strcmp(atype_string, "major") == 0) {
        *atype = ALLELE_MAJOR;
    }
    else if (strcmp(atype_string, "nonmajor") == 0) {
        *atype = ALLELE_NONMAJOR;
    }
    else {
        error("Error: allele type (%s) not recognised. Must be one of nref|alt1|minor|major|nonmajor: %s\n", atype_string);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools view [options] <in.vcf.gz> [region1 [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -G,   --drop-genotypes              drop individual genotype information (after subsetting if -s option set)\n");
    fprintf(stderr, "    -h/H, --header-only/--no-header     print the header only/suppress the header in VCF output\n");
    fprintf(stderr, "    -l,   --compression-level [0-9]     compression level: 0 uncompressed, 1 best speed, 9 best compression [%d]\n", args->clevel);
    fprintf(stderr, "    -o,   --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>              restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>           restrict to regions listed in a file\n");
    fprintf(stderr, "    -t, --targets [^]<region>           similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n");
    fprintf(stderr, "    -T, --targets-file [^]<file>        similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Subset options:\n");
    fprintf(stderr, "    -a, --trim-alt-alleles        trim alternate alleles not seen in the subset\n");
    fprintf(stderr, "    -I, --no-update               do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)\n");
    fprintf(stderr, "    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "    -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "        --force-samples           only warn about unknown subset samples\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -c/C, --min-ac/--max-ac <int>[:<type>]      minimum/maximum count for non-reference (nref), 1st alternate (alt1), least frequent\n");
    fprintf(stderr, "                                                   (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n");
    fprintf(stderr, "    -f,   --apply-filters <list>                require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes\n");
    fprintf(stderr, "    -i/e, --include/--exclude <expr>            select/exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -k/n, --known/--novel                       select known/novel sites only (ID is not/is '.')\n");
    fprintf(stderr, "    -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)\n");
    fprintf(stderr, "    -p/P, --phased/--exclude-phased             select/exclude sites where all samples are phased\n");
    fprintf(stderr, "    -q/Q, --min-af/--max-af <float>[:<type>]    minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent\n");
    fprintf(stderr, "                                                   (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n");
    fprintf(stderr, "    -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype\n");
    fprintf(stderr, "    -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,other [null]\n");
    fprintf(stderr, "    -x/X, --private/--exclude-private           select/exclude sites where the non-reference alleles are exclusive (private) to the subset samples\n");
    fprintf(stderr, "\n");
    exit(1);
}

#define ZERO_BASED_POSITION 1
#define ZERO_BASED_SAMPLE_IDS 1
unsigned char g_CSV_MISSING_CHARACTER='#';
#define CSV_NON_REF_REPRESENTATION "&"

unsigned string_to_tiledb_output_version(const char* string)
{
    if(strcmp(string, "v0") == 0)
        return TILEDB_OUTPUT_CSV_V0;
    else
        if(strcmp(string, "csv_v1") == 0)
            return  TILEDB_OUTPUT_CSV_V1;
        else
            if(strcmp(string, "bin_v1") == 0)
                return TILEDB_OUTPUT_BINARY_V1;
            else
                error("Unknown TileDB output format: '%s'\n",string);
    exit(-1);
    return 0;
}

void tiledb_csv_printer(buffer_wrapper* buffer, unsigned bcf_ht_type, const char* fmt_string, ...)
{
    va_list ap;
    va_start(ap, fmt_string);
    int buffer_space_available = buffer->m_size - buffer->m_offset;
    //If fmt_string is NULL, print the MISSING value
    int num_chars_printed = fmt_string ? vsnprintf(buffer->m_buffer+buffer->m_offset, buffer_space_available, fmt_string, ap)
        : snprintf(buffer->m_buffer+buffer->m_offset, buffer_space_available, ",%c", g_CSV_MISSING_CHARACTER);
    if(num_chars_printed >= buffer_space_available) 
    { 
        fprintf(stderr, "TileDB CSV line too long - exiting :(\n"); 
        exit(-1); 
    } 
    buffer->m_offset += num_chars_printed; 
    va_end(ap);
}

void tiledb_binary_printer(buffer_wrapper* buffer, unsigned bcf_ht_type, const char* fmt_string, ...)
{
    va_list ap;
    va_start(ap, fmt_string);
    //Temporary variables for use in switch statement
    int int_val;
    float float_val;
    char char_val;
    int64_t int64_val;
    //To be used in memcpy
    size_t bytes = 0;
    char* ptr;
    char missing_string[] = { g_CSV_MISSING_CHARACTER, '\0' };
    //For proper printing
    //If fmt_string is NULL, print missing value
#define TILEDB_BIN_SET_PTR(va_type_t, type_t, val, missing_value)\
    val = fmt_string ? (type_t)(va_arg(ap, va_type_t)) : missing_value; \
    ptr = (char*)(&(val)); \
    bytes = sizeof(type_t);
    switch(bcf_ht_type)
    {
        case BCF_HT_INT:
            TILEDB_BIN_SET_PTR(int, int, int_val, INT_MAX);
            break;
        case BCF_HT_REAL:
            TILEDB_BIN_SET_PTR(double, float, float_val, FLT_MAX);
            break;
        case BCF_HT_STR:
            ptr = fmt_string ? va_arg(ap, char*) : missing_string;
            bytes = strlen(ptr);
            break;
        case BCF_HT_CHAR:
            TILEDB_BIN_SET_PTR(int, char, char_val, g_CSV_MISSING_CHARACTER);
            break;
        case BCF_HT_INT64:
            TILEDB_BIN_SET_PTR(int64_t, int64_t, int64_val, LLONG_MAX);
            break;
        case BCF_HT_VOID:
            bytes = 0;
            break;
        default:
            fprintf(stderr, "Unknown BCF HT type %d, exiting\n", bcf_ht_type);
            exit(-1);
    }
    va_end(ap);
    if(bytes == 0)
        return;
    int buffer_space_available = buffer->m_size - buffer->m_offset;
    if(bytes > buffer_space_available) 
    { 
        fprintf(stderr, "TileDB binary line too long - exiting :(\n"); 
        exit(-1); 
    }
    memcpy(buffer->m_buffer + buffer->m_offset, ptr, bytes); 
    buffer->m_offset += bytes; 
}

typedef struct
{
  int input_idx;
  int64_t contig_length;
  int64_t* input_2_global_idx_array;
  int64_t* contig_offsets;
}sqlite_data_struct;

int sqlite_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 1);
  int64_t global_idx = strtoll(field_values[0], 0, 10);
  sqlite_data_struct* data = (sqlite_data_struct*)ptr;
  data->input_2_global_idx_array[data->input_idx] = global_idx;
  return 0;
}

int contig_sqlite_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
    assert(num_columns == 3);
    if(sqlite_handler(ptr, 1, field_values, column_names) != 0)
        return -1;
    sqlite_data_struct* data = (sqlite_data_struct*)ptr;
    data->contig_length = strtoull(field_values[1], 0, 10);
    data->contig_offsets[data->input_idx] = strtoull(field_values[2], 0, 10);
    return 0;
}

int sqlite_name_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 1);
  char* name = (char*)ptr;
  strcpy(name, field_values[0]);
  return 0;
}

typedef struct
{
    int64_t contig_offset;
    int64_t contig_length;
}max_contig_data_struct;

int max_contig_offset_sqlite_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
    assert(num_columns == 2);
    max_contig_data_struct* data = (max_contig_data_struct*)ptr;
    data->contig_length = strtoll(field_values[0], 0, 10);
    data->contig_offset = strtoll(field_values[1], 0, 10);
    return 0;
}

#define COMMON_QUERY_VARIABLE_DECLARATION \
    sqlite_mappings_struct* mapping_info = (sqlite_mappings_struct*)info_ptr; \
    char* error_msg = 0; \
    char query_string[4096]; \
    char insert_string[4096]; \
    sqlite_data_struct data; \
    int i = 0;

void query_sample_name(void* info_ptr, int64_t sample_idx, char* sample_name)
{
    sqlite_mappings_struct* mapping_info = (sqlite_mappings_struct*)info_ptr;
    char* error_msg = 0;
    char query_string[4096];
#ifdef ZERO_BASED_SAMPLE_IDS
    ++sample_idx;
#endif
    sample_name[0] = '\0';
    sprintf(query_string,"select sample_name from sample_names where sample_names.sample_idx == %"PRIi64";",sample_idx);
    sqlite3_exec(mapping_info->db, query_string, sqlite_name_handler, sample_name, &error_msg);
    assert(sample_name[0] != '\0');
}

typedef struct
{
    int64_t m_contig_offset;
    char* m_contig_name;
} contig_name_offset_struct;

int sqlite_contig_name_offset_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 2);
  contig_name_offset_struct* data = (contig_name_offset_struct*)ptr;
  strcpy(data->m_contig_name, field_values[0]);
  data->m_contig_offset = strtoll(field_values[1], 0, 10);
  return 0;
}

int64_t query_contig_name(void* info_ptr, int64_t column_idx, char* contig_name)
{
    sqlite_mappings_struct* mapping_info = (sqlite_mappings_struct*)info_ptr;
    char* error_msg = 0;
    char query_string[4096];
    contig_name_offset_struct data;
    data.m_contig_name = contig_name;
    data.m_contig_offset = column_idx;
#ifndef ZERO_BASED_POSITION
    --column_idx;
#endif
    contig_name[0] = '\0';
    sprintf(query_string,"select contig_name,contig_offset from contig_names where contig_names.contig_offset<=%"PRIi64" and contig_names.contig_offset+contig_names.contig_length > %"PRIi64";",column_idx, column_idx);
    sqlite3_exec(mapping_info->db, query_string, sqlite_contig_name_offset_handler, &data, &error_msg);
    assert(contig_name[0] != '\0');
    return (column_idx - data.m_contig_offset);
}

const int64_t* query_samples_idx(void* info_ptr, int n_samples, const char* const* sample_names)
{
    COMMON_QUERY_VARIABLE_DECLARATION; 
    //Memory allocation for mappings
    if(n_samples)
    {
        mapping_info->input_sample_idx_2_global_idx = (int64_t*)realloc(mapping_info->input_sample_idx_2_global_idx, n_samples*sizeof(int64_t));
        memset(mapping_info->input_sample_idx_2_global_idx, -1, n_samples*sizeof(int64_t));  //initialize to -1
    }
    //Samples 
    data.input_2_global_idx_array = mapping_info->input_sample_idx_2_global_idx;
    for(i=0;i<n_samples;++i)
    {
        data.input_idx = i; 
        sprintf(query_string,"select sample_idx from sample_names where sample_names.sample_name == \"%s\";",sample_names[i]);
        sqlite3_exec(mapping_info->db, query_string, sqlite_handler, &data, &error_msg);
        if(data.input_2_global_idx_array[i] == -1)
        {
            sprintf(insert_string,"insert into sample_names (sample_name) values( \"%s\" );",sample_names[i]);
            sqlite3_exec(mapping_info->db, insert_string, 0, 0, &error_msg);
            sqlite3_exec(mapping_info->db, query_string, sqlite_handler, &data, &error_msg);
            assert(data.input_2_global_idx_array[i] >= 0);
        }
#ifdef ZERO_BASED_SAMPLE_IDS
        --(mapping_info->input_sample_idx_2_global_idx[i]);
#endif
    } 
    return mapping_info->input_sample_idx_2_global_idx;
}

const int64_t* query_contigs_offset(void* info_ptr, int n_contigs, const char* const* contig_names, const int64_t* contig_lengths)
{
    COMMON_QUERY_VARIABLE_DECLARATION; 
    //max contig offset
    const char* max_contig_offset_query = "select contig_length, contig_offset from contig_names order by contig_offset desc limit 1;";
    max_contig_data_struct max_contig;
    //First time
    max_contig.contig_offset = 0ull;
    max_contig.contig_length = 0ull;
    //Memory allocation for mappings
    if(n_contigs)
    {
        mapping_info->input_contig_idx_2_global_idx = (int64_t*)realloc(mapping_info->input_contig_idx_2_global_idx, n_contigs*sizeof(int64_t));
        memset(mapping_info->input_contig_idx_2_global_idx, -1, n_contigs*sizeof(int64_t));     //initialize to -1
        mapping_info->input_contig_idx_2_offset = (int64_t*)realloc(mapping_info->input_contig_idx_2_offset, n_contigs*sizeof(uint64_t));
    }
    //Contigs
    data.input_2_global_idx_array = mapping_info->input_contig_idx_2_global_idx;
    data.contig_offsets = mapping_info->input_contig_idx_2_offset;
    for(i=0;i<n_contigs;++i)
    {
        //query
        data.input_idx = i; 
        sprintf(query_string,"select contig_idx,contig_length,contig_offset from contig_names where contig_names.contig_name == \"%s\";",
                contig_names[i]);
        sqlite3_exec(mapping_info->db, query_string, contig_sqlite_handler, &data, &error_msg);
        //Contig does not exist in DB
        if(data.input_2_global_idx_array[i] == -1)
        {
            assert(contig_lengths);     //should be non-NULL
            uint64_t curr_contig_length = contig_lengths[i];
            assert(curr_contig_length > 0);
            //why the retry?
            //the schema forces contig_offset to have unique values, if multiple inserts are tried in parallel, only one will succeed 
            while(data.input_2_global_idx_array[i] == -1)
            {
                //get max offset value from SQL
                sqlite3_exec(mapping_info->db, max_contig_offset_query, max_contig_offset_sqlite_handler, &max_contig, &error_msg);
                //Compute offset for new entry and insert
                sprintf(insert_string,"insert into contig_names (contig_name, contig_length, contig_offset) values( \"%s\",%" PRIu64 ",%" PRIu64 " );",
                        contig_names[i], curr_contig_length, max_contig.contig_offset + max_contig.contig_length);
                error_msg = 0;
                //May fail because of parallel insert
                if(sqlite3_exec(mapping_info->db, insert_string, 0, 0, &error_msg) == SQLITE_OK)
                {
                    sqlite3_exec(mapping_info->db, query_string, contig_sqlite_handler, &data, &error_msg);
                    assert(data.input_2_global_idx_array[i] >= 0);
                }
                else
                    if(error_msg)
                        free(error_msg);
            }
        }
        if(contig_lengths)
            assert(data.contig_length == contig_lengths[i]);       //should be same as current contig
    }
    return mapping_info->input_contig_idx_2_offset;
}

const int64_t* query_fields_idx(void* info_ptr, int n_fields, const char* const* field_names)
{
    COMMON_QUERY_VARIABLE_DECLARATION; 
    if(n_fields)
    {
        mapping_info->input_field_idx_2_global_idx = (int64_t*)realloc(mapping_info->input_field_idx_2_global_idx, n_fields*sizeof(int64_t));
        memset(mapping_info->input_field_idx_2_global_idx, -1,n_fields*sizeof(int64_t));        //initialize to -1 
    }
    //INFO/FILTER/FORMAT fields
    data.input_2_global_idx_array = mapping_info->input_field_idx_2_global_idx;
    for(i=0;i<n_fields;++i)
    {
        data.input_idx = i; 
        sprintf(query_string,"select field_idx from field_names where field_names.field_name == \"%s\";", field_names[i]);
        sqlite3_exec(mapping_info->db, query_string, sqlite_handler, &data, &error_msg);
        if(data.input_2_global_idx_array[i] == -1)
        {
            sprintf(insert_string,"insert into field_names (field_name) values( \"%s\" );", field_names[i]);
            sqlite3_exec(mapping_info->db, insert_string, 0, 0, &error_msg);
            sqlite3_exec(mapping_info->db, query_string, sqlite_handler, &data, &error_msg);
            assert(data.input_2_global_idx_array[i] >= 0);
        }
    }
    return mapping_info->input_field_idx_2_global_idx;
}

void initialize_samples_contigs_and_fields_idx(sqlite_mappings_struct* mapping_info, const bcf_hdr_t* hdr)
{
    int n_samples = bcf_hdr_nsamples(hdr);
    int n_fields = hdr->n[BCF_DT_ID];
    int n_contigs = hdr->n[BCF_DT_CTG];
    char** ptr = 0;
    //Allocate arrays to store field and contig names
    mapping_info->m_field_names = (char**)malloc(n_fields*sizeof(char*));
    mapping_info->m_contig_names = (char**)malloc(n_contigs*sizeof(char*));
    mapping_info->m_contig_lengths = (int64_t*)malloc(n_contigs*sizeof(uint64_t));
    //Initialize length variables
    mapping_info->m_num_samples = n_samples;
    mapping_info->m_num_fields = n_fields;
    mapping_info->m_num_contigs = n_contigs;

    int i = 0;
    //Contigs
    ptr = mapping_info->m_contig_names;
    for(i=0;i<hdr->n[BCF_DT_CTG];++i)
    {
        ptr[i] = (char*)(bcf_hdr_int2id(hdr, BCF_DT_CTG, i));
        mapping_info->m_contig_lengths[i] = bcf_hdr_id2contig_length(hdr, i); 
    }
    //Fields
    ptr = mapping_info->m_field_names;
    for(i=0;i<hdr->n[BCF_DT_ID];++i)
        ptr[i] = (char*)(bcf_hdr_int2id(hdr, BCF_DT_ID, i));
    //Query samples, contigs and fields
    if(mapping_info->m_tiledb_override_sample_name)
    {
        assert(n_samples == 1); //header should have one sample, if the sample name is overridden by a command line argument
        query_samples_idx(mapping_info, 1, (const char* const*)(&(mapping_info->m_tiledb_override_sample_name)));
    }
    else
        query_samples_idx(mapping_info, n_samples, (const char* const*)(hdr->samples));
    query_contigs_offset(mapping_info, n_contigs, (const char* const*)(mapping_info->m_contig_names), mapping_info->m_contig_lengths);
    query_fields_idx(mapping_info, n_fields, (const char* const*)(mapping_info->m_field_names));
    //Version dependent values
    switch(mapping_info->m_tiledb_output_version)
    {
        case TILEDB_OUTPUT_CSV_V0:
            g_CSV_MISSING_CHARACTER = '#';
            mapping_info->m_tiledb_printer = tiledb_csv_printer; 
            break;
        case TILEDB_OUTPUT_CSV_V1:
            g_CSV_MISSING_CHARACTER = '*';
            mapping_info->m_tiledb_printer = tiledb_csv_printer; 
            break;
        case TILEDB_OUTPUT_BINARY_V1:
            g_CSV_MISSING_CHARACTER = '\0';
            mapping_info->m_tiledb_printer = tiledb_binary_printer; 
            break;
        default:
            error("Unknown TileDB output format %d\n", mapping_info->m_tiledb_output_version);
            break;
    }
}

int sqlite_count_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 1);
  //Ptr to count variable
  int64_t* count_ptr = (int64_t*)ptr;
  *count_ptr = strtoll(field_values[0], 0, 10);
  return 0;
}

int sqlite_scan_samples_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 2);
  sqlite_mappings_struct* info = (sqlite_mappings_struct*)ptr;
  int64_t sample_idx = strtoll(field_values[0], 0, 10);
#ifdef ZERO_BASED_SAMPLE_IDS
  //SQLite idxs are 1 based, but zero based in TileDB
  --sample_idx;
#endif
  assert(sample_idx < info->m_num_samples);
  info->m_sample_names[sample_idx] = strdup(field_values[1]);
  return 0;
}

void read_all_samples_from_sqlite(sqlite_mappings_struct* mapping_info)
{
    char* error_msg = 0;
    //Find #samples in SQLite
    const char* query_string = "select count(*) from sample_names";
    sqlite3_exec(mapping_info->db, query_string, sqlite_count_handler, &(mapping_info->m_num_samples), &error_msg);
#ifndef ZERO_BASED_SAMPLE_IDS
    //SQLite idxs are 1 based, so increase size of array if NOT 0-based
    ++(mapping_info->m_num_samples);
#endif
    //allocate array
    mapping_info->m_sample_names = (char**)calloc(mapping_info->m_num_samples, sizeof(char*));
    //scan table
    const char* scan_query_string = "select sample_idx,sample_name from sample_names";
    sqlite3_exec(mapping_info->db, scan_query_string, sqlite_scan_samples_handler, mapping_info, &error_msg);
    //flag that all strings are allocated and must be freed later
    mapping_info->m_sample_name_strings_allocated = 1;
}

int sqlite_scan_fields_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 2);
  sqlite_mappings_struct* info = (sqlite_mappings_struct*)ptr;
  int64_t field_idx = strtoll(field_values[0], 0, 10);
  assert(field_idx < info->m_num_fields);
  info->m_field_names[field_idx] = strdup(field_values[1]);
  return 0;
}

void read_all_fields_from_sqlite(sqlite_mappings_struct* mapping_info)
{
    char* error_msg = 0;
    //Find #fields in SQLite
    const char* query_string = "select count(*) from field_names";
    sqlite3_exec(mapping_info->db, query_string, sqlite_count_handler, &(mapping_info->m_num_fields), &error_msg);
    //SQLite idxs are 1 based, so increase size of array
    ++(mapping_info->m_num_fields);
    //allocate array
    mapping_info->m_field_names = (char**)calloc(mapping_info->m_num_fields, sizeof(char*));
    //scan table
    const char* scan_query_string = "select field_idx,field_name from field_names";
    sqlite3_exec(mapping_info->db, scan_query_string, sqlite_scan_fields_handler, mapping_info, &error_msg);
    //flag that all strings are allocated and must be freed later
    mapping_info->m_field_name_strings_allocated = 1;
}

int sqlite_scan_contigs_handler(void* ptr, int num_columns, char** field_values, char** column_names)
{
  assert(num_columns == 4);
  sqlite_mappings_struct* info = (sqlite_mappings_struct*)ptr;
  int64_t contig_idx = strtoll(field_values[0], 0, 10);
  assert(contig_idx < info->m_num_contigs);
  info->m_contig_names[contig_idx] = strdup(field_values[1]);
  info->input_contig_idx_2_offset[contig_idx] = strtoll(field_values[2], 0, 10);
  info->m_contig_lengths[contig_idx] = strtoll(field_values[3], 0, 10);
  return 0;
}

void read_all_contigs_from_sqlite(sqlite_mappings_struct* mapping_info)
{
    char* error_msg = 0;
    //Find #contigs in SQLite
    const char* query_string = "select count(*) from contig_names";
    sqlite3_exec(mapping_info->db, query_string, sqlite_count_handler, &(mapping_info->m_num_contigs), &error_msg);
    //SQLite idxs are 1 based, so increase size of array
    ++(mapping_info->m_num_contigs);
    //allocate arrays
    mapping_info->m_contig_names = (char**)calloc(mapping_info->m_num_contigs, sizeof(char*));
    mapping_info->m_contig_lengths = (int64_t*)calloc(mapping_info->m_num_contigs, sizeof(int64_t));
    //initialize to -1, invalid
    memset(mapping_info->m_contig_lengths, -1, mapping_info->m_num_contigs*sizeof(int64_t));
    mapping_info->input_contig_idx_2_offset = (int64_t*)calloc(mapping_info->m_num_contigs, sizeof(int64_t));
    //initialize to -1, invalid
    memset(mapping_info->input_contig_idx_2_offset, -1, mapping_info->m_num_contigs*sizeof(int64_t));
    //scan table
    const char* scan_query_string = "select contig_idx,contig_name,contig_offset,contig_length from contig_names";
    sqlite3_exec(mapping_info->db, scan_query_string, sqlite_scan_contigs_handler, mapping_info, &error_msg);
    //flag that all strings are allocated and must be freed later
    mapping_info->m_contig_name_strings_allocated = 1;
}

//Reads all tables from SQLite into memory
void read_all_from_sqlite(sqlite_mappings_struct* mapping_info)
{
    read_all_samples_from_sqlite(mapping_info);
    read_all_fields_from_sqlite(mapping_info);
    read_all_contigs_from_sqlite(mapping_info); 
}

void test_read_all(const char* sqlite_file)
{
    sqlite_mappings_struct* info = (sqlite_mappings_struct*)allocate_sqlite3_mapping(sqlite_file);
    read_all_from_sqlite(info);
    int64_t i = 0;
    for(i=0;i<info->m_num_samples;++i)
        printf("Idx %"PRIi64" sample name %s\n", i, info->m_sample_names[i]);
    for(i=0;i<info->m_num_fields;++i)
        printf("Idx %"PRIi64" field name %s\n", i, info->m_field_names[i]);
    for(i=0;i<info->m_num_contigs;++i)
        printf("Idx %"PRIi64" contig name %s offset %"PRIi64" length %"PRIi64"\n", i, info->m_contig_names[i],
                info->input_contig_idx_2_offset[i], info->m_contig_lengths[i]);
    free_sqlite3_data(info);
    free(info);
}

int write_csv_line(sqlite_mappings_struct* mapping_info, csv_output_struct* csv_output_info,
	bcf_hdr_t* out_hdr, bcf1_t* line, int input_sample_idx)
{
    //Buffer to write CSV line
    buffer_wrapper* csv_out_buffer = &(csv_output_info->m_csv_out_buffer);
    //Sampling info
    random_sampling_struct* sampling_info = SHOULD_DO_SAMPLING(csv_output_info->m_sampling_info) ? &(csv_output_info->m_sampling_info) : 0 ;
    //Profiling info
    gvcf_stat_struct* profile_intervals = SHOULD_DO_PROFILING(csv_output_info->m_profile_intervals) ? &(csv_output_info->m_profile_intervals) : 0; 
    //Buffer to use to call get_info/format_values()
    char** buffer = &(csv_output_info->m_htslib_buffer);
    unsigned* buffer_size = &(csv_output_info->m_htslib_buffer_size);
    //Printer function
    TileDBPrinterTy TILEDB_CSV_BPRINTF = mapping_info->m_tiledb_printer;
    int contig_id = line->rid;
    ASSERT(mapping_info->input_contig_idx_2_global_idx[contig_id] >= 0);
    int64_t contig_offset = mapping_info->input_contig_idx_2_offset[contig_id];
    //old version
    uint8_t is_old_tiledb_version = (mapping_info->m_tiledb_output_version == TILEDB_OUTPUT_CSV_V0);
    uint8_t is_tiledb_binary_output = (mapping_info->m_tiledb_output_version == TILEDB_OUTPUT_BINARY_V1);
    uint8_t is_anchor_cell = csv_output_info->m_print_tiledb_anchor_cell;
    unsigned offset_at_start = csv_out_buffer->m_offset;
    unsigned cell_size_offset = 0u;
    {
        int num_values = -1;
#define PRINT_TILEDB_CSV_WITH_GT_FLAG(field_name, bcf_ht_type, get_function, type_t, format_specifier, vector_end_condition, missing_condition, length_descriptor, fixed_num_values, is_gt_field) \
        { \
            int num_elements = (*buffer_size)/sizeof(type_t); \
            /*For anchor cells, print NULL only*/ \
            num_values = is_anchor_cell ? -1 : get_function(out_hdr, line, field_name, (void**)buffer, &num_elements, bcf_ht_type); \
            if(num_elements*sizeof(type_t) > (*buffer_size)) \
                (*buffer_size) = num_elements*sizeof(type_t); \
            if(num_values < 0) \
            { \
                /*For newer versions, for attributes with variable# elements, print #elements first*/ \
                if(length_descriptor != BCF_VL_FIXED && !is_old_tiledb_version) \
                    TILEDB_CSV_BPRINTF(csv_out_buffer,bcf_ht_type,",0"); \
                else \
                { \
                    int i=0; \
                    for(i=0;i<(fixed_num_values);++i) \
                        TILEDB_CSV_BPRINTF(csv_out_buffer,bcf_ht_type,0,g_CSV_MISSING_CHARACTER);  \
                } \
            } \
            else \
            { \
                /*If GT field is being printed or (for newer versions and for attributes with variable# elements) print #elements first*/ \
                if(is_gt_field || (length_descriptor != BCF_VL_FIXED && !is_old_tiledb_version)) \
                    TILEDB_CSV_BPRINTF(csv_out_buffer,bcf_ht_type,",%d",num_values); \
                int k = 0;  \
                type_t* ptr = (type_t*)(*buffer);      \
                for(k=0;k<num_values;++k)   \
                {   \
                    type_t val = ptr[k];    \
                    if(is_gt_field) \
                        val = bcf_gt_allele((int)val); \
                    TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_VOID,","); \
                    if(missing_condition) \
                        TILEDB_CSV_BPRINTF(csv_out_buffer,bcf_ht_type,0,g_CSV_MISSING_CHARACTER);  \
                    TILEDB_CSV_BPRINTF(csv_out_buffer,bcf_ht_type,format_specifier,val);  \
                }   \
            } \
        }

        //This macro does not have a is_gt_field flag, is_gt_field assumed to be false
#define PRINT_TILEDB_CSV(field_name, bcf_ht_type, get_function, type_t, format_specifier, vector_end_condition, missing_condition, length_descriptor, fixed_num_values) \
        PRINT_TILEDB_CSV_WITH_GT_FLAG(field_name, bcf_ht_type, get_function, type_t, format_specifier, vector_end_condition, missing_condition, length_descriptor, fixed_num_values, 0)

        //Print co-ordinates : sample id, location
        ASSERT(mapping_info->input_sample_idx_2_global_idx[input_sample_idx] >= 0);
#ifdef ZERO_BASED_POSITION
        int64_t position = contig_offset + ((uint64_t)line->pos);
#else
        int64_t position = contig_offset + ((uint64_t)line->pos) + 1;
#endif
        //For anchor cells, the position is the anchor position
        if(is_anchor_cell)
            position = csv_output_info->m_curr_tiledb_anchor_position;
	if(sampling_info)
	{
	    int random_value = rand();
	    if(random_value < sampling_info->m_sampling_limit)
		fprintf(sampling_info->m_spit_random_positions,"%"PRIi64"\n",position);
	    return 0;
	}
        //get END position
	int64_t global_end_position = 0;
        if(is_anchor_cell)
            global_end_position = -1;   //why -1, essentially this tells TileDB that this is an invalid cell (anchor cell)
        else
        {
            int num_elements = (*buffer_size)/sizeof(int);
            int num_values = bcf_get_info_values(out_hdr, line, "END", (void**)buffer, &num_elements, BCF_HT_INT);
            //END not found, single position record
            if(num_values < 0)
            {
                global_end_position = position;
                if(csv_output_info->m_treat_deletions_as_intervals)
                {
                    //Check for deletions
                    int j = 0;
                    char** alleles = line->d.allele;
                    int ref_length = strlen(alleles[0]);
                    for(j=1;j<line->n_allele;++j)
                    {
                        if(bcf_get_variant_type(line, j) == VCF_INDEL
                                && ref_length > strlen(alleles[j]))
                        {
                            global_end_position = position+ref_length-1;
                            break;
                        }
                    }
                }
            }
            else
            {
                int end_pos = ((int*)(*buffer))[0];
#ifdef ZERO_BASED_POSITION
                --end_pos;
#endif
		global_end_position = contig_offset + ((int64_t)end_pos);
            }
            csv_output_info->m_last_inserted_tiledb_position = global_end_position;
        }
	if(profile_intervals)
	{
	    uint64_t interval_length = (global_end_position - position + 1);
	    profile_intervals->m_num_valid_positions += interval_length;
	    profile_intervals->m_sum_sq_valid_interval_length += (interval_length*interval_length);
	    if(interval_length > profile_intervals->m_max_valid_interval_length)
		profile_intervals->m_max_valid_interval_length = interval_length;
	    ++(profile_intervals->m_num_valid_intervals);
	    uint64_t last_invalid_interval_length = position - profile_intervals->m_last_valid_interval_end - 1;
	    profile_intervals->m_last_valid_interval_end = global_end_position;
	    if(last_invalid_interval_length > 0)
	    {
		profile_intervals->m_num_invalid_positions += last_invalid_interval_length;
		profile_intervals->m_sum_sq_invalid_interval_length += (last_invalid_interval_length*last_invalid_interval_length);
		if(last_invalid_interval_length > profile_intervals->m_max_invalid_interval_length)
		    profile_intervals->m_max_invalid_interval_length = last_invalid_interval_length;
		++(profile_intervals->m_num_invalid_intervals);
	    }
	    if(profile_intervals->m_max_position < position)
		profile_intervals->m_max_position = position;
	    return 0;
	}
        csv_out_buffer->m_global_sample_idx = mapping_info->input_sample_idx_2_global_idx[input_sample_idx];
        csv_out_buffer->m_global_start_position = position;
        csv_out_buffer->m_global_end_position = global_end_position;
        if(!(csv_output_info->m_skip_coordinates))
        {
            //print sample idx and position
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT64,"%"PRIi64, mapping_info->input_sample_idx_2_global_idx[input_sample_idx]);
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT64,",%"PRIi64, position);
            //For binary output, size of total cell is required after co-ordinates
            //FIXME:The type of size is size_t, assuming 64-bit on 64-bit Linux platforms
            //Just allocate space for now and update later once the cell size is known
            if(is_tiledb_binary_output)
            {
                cell_size_offset = csv_out_buffer->m_offset;
                ALLOCATE_SPACE_IN_BUFFER(csv_out_buffer, sizeof(int64_t));
            }
            //print END position
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT64,",%"PRIi64, global_end_position);
        }
        //For anchor cells, print NULL for REF, ALT, QUAL, FILTER fields
        if(is_anchor_cell)
        {
            //REF
            //For binary, the size of the var:char field is all that is needed 
            if(is_tiledb_binary_output)
                TILEDB_CSV_BPRINTF(csv_out_buffer, BCF_HT_INT,",%d", 0);
            else
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_STR,",%c", g_CSV_MISSING_CHARACTER);
            //ALT
            //For binary, the size of the var:char field is all that is needed 
            if(is_tiledb_binary_output)
                TILEDB_CSV_BPRINTF(csv_out_buffer, BCF_HT_INT,",%d", 0);
            else
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_STR,",%c", g_CSV_MISSING_CHARACTER);
            //QUAL
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_REAL,0,g_CSV_MISSING_CHARACTER);
            //Set #filters to 0
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT,",%d",0);
        }
        else
        {
            //For binary, the size of the var:char field must be printed first 
            //Print length of REF string 
            if(is_tiledb_binary_output)
                TILEDB_CSV_BPRINTF(csv_out_buffer, BCF_HT_INT,",%d", strlen(line->d.allele[0]));
            //print reference allele
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_STR,",%s", line->d.allele[0]);
            //Position in buffer after REF
            int after_ref_offset = csv_out_buffer->m_offset;
            //Separator for ALT alleles
            char alt_separator = ',';
            if(is_old_tiledb_version)
            {
                //For old versions, print #ALT alleles first
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT,",%d", line->n_allele-1);
                alt_separator = ',';
            }
            else
                alt_separator = '|';
            //For binary, the size of the var:char field must be printed first 
            //Just allocate space for now and update later once the length is known
            if(is_tiledb_binary_output)
            {
                ALLOCATE_SPACE_IN_BUFFER(csv_out_buffer, sizeof(int));
            }
            int alt_start_offset = csv_out_buffer->m_offset;
            //print alt alleles
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_VOID,",");   //print comma before first ALT allele
            int j = 0;
            for(j=1;j<line->n_allele;++j)
            {
                if(j > 1)
                    TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_CHAR,"%c",alt_separator);
                if(strcmp(line->d.allele[j],"<NON_REF>") == 0)
                    TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_STR,"%s",CSV_NON_REF_REPRESENTATION);
                else
                    TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_STR,"%s",line->d.allele[j]);
            }
            //For binary output, print length of concatenated ALT string 
            if(is_tiledb_binary_output)
            {
                int alt_length = csv_out_buffer->m_offset - alt_start_offset;
                memcpy(csv_out_buffer->m_buffer+after_ref_offset, &alt_length, sizeof(int));
            }
            //print QUAL
            if(bcf_float_is_missing(line->qual))
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_REAL,0,g_CSV_MISSING_CHARACTER);
            else
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_REAL,",%0.3f",line->qual);
            //print number of filters, followed by filter idx
            TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT,",%d",line->d.n_flt);
            for(j=0;j<line->d.n_flt;++j)
            {
                ASSERT(line->d.flt[j] >= 0 && line->d.flt[j] < out_hdr->n[BCF_DT_ID]);
                ASSERT(mapping_info->input_field_idx_2_global_idx[line->d.flt[j]] >= 0);
                TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_INT64,",%"PRIi64,mapping_info->input_field_idx_2_global_idx[line->d.flt[j]]);
            }
        }
        //Print relevant INFO fields
        PRINT_TILEDB_CSV("BaseQRankSum",BCF_HT_REAL, bcf_get_info_values, float,"%f",0,(bcf_float_is_missing(val)),BCF_VL_FIXED,1);
        PRINT_TILEDB_CSV("ClippingRankSum",BCF_HT_REAL, bcf_get_info_values,float,"%f",0,(bcf_float_is_missing(val)),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("MQRankSum",BCF_HT_REAL, bcf_get_info_values,float,"%f",0,(bcf_float_is_missing(val)),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("ReadPosRankSum",BCF_HT_REAL, bcf_get_info_values,float,"%f",0,(bcf_float_is_missing(val)),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("DP",BCF_HT_INT, bcf_get_info_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("MQ",BCF_HT_REAL, bcf_get_info_values,float,"%f",0,(bcf_float_is_missing(val)),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("MQ0",BCF_HT_INT, bcf_get_info_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 1);
        //Print relevant FORMAT fields
        PRINT_TILEDB_CSV("DP",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("MIN_DP",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("GQ",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 1);
        PRINT_TILEDB_CSV("SB",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_FIXED, 4);
        PRINT_TILEDB_CSV("AD",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_R, line->n_allele);
        int num_gts =(line->n_allele*(line->n_allele+1))/2;
        PRINT_TILEDB_CSV("PL",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_G, num_gts);
        PRINT_TILEDB_CSV_WITH_GT_FLAG("GT",BCF_HT_INT, bcf_get_format_values,int32_t,"%d",0,(val == bcf_int32_missing),BCF_VL_P, 0, 1);
    }
    TILEDB_CSV_BPRINTF(csv_out_buffer,BCF_HT_VOID,"\n");
    //For binary output, print cell size at the correct spot
    //FIXME: assuming 64 bit size_t on Linux 64 bit
    if(is_tiledb_binary_output && !(csv_output_info->m_skip_coordinates))
    {
        int64_t cell_size = csv_out_buffer->m_offset - offset_at_start;
        memcpy(csv_out_buffer->m_buffer+cell_size_offset, &cell_size, sizeof(cell_size));
    }
    unsigned curr_offset = csv_out_buffer->m_offset; 
    //Check if need to write to stream
    FILE* csv_fptr = csv_output_info->m_csv_out_fptr;
    if(csv_fptr)
    {
      fputs(csv_out_buffer->m_buffer, csv_fptr);
      csv_out_buffer->m_offset = 0;
    }
    return curr_offset;
}

void write_csv_for_one_VCF_line(sqlite_mappings_struct* mapping_info, csv_output_struct* csv_output_info,
	bcf_hdr_t* out_hdr, bcf1_t* line)
{
    int i = 0;
    bcf_unpack(line, BCF_UN_ALL);
    if(SHOULD_DO_SAMPLING(csv_output_info->m_sampling_info) || SHOULD_DO_PROFILING(csv_output_info->m_profile_intervals))
	write_csv_line(mapping_info, csv_output_info, out_hdr, line, 0);
    else
    {
        //Insert anchor cells if needed
        if(csv_output_info->m_tiledb_anchor_interval)
        {
            int contig_id = line->rid;
            ASSERT(mapping_info->input_contig_idx_2_global_idx[contig_id] >= 0);
            int64_t contig_offset = mapping_info->input_contig_idx_2_offset[contig_id];
#ifdef ZERO_BASED_POSITION
            int64_t position = contig_offset + ((uint64_t)line->pos);
#else
            int64_t position = contig_offset + ((uint64_t)line->pos) + 1;
#endif
            int64_t curr_anchor_position = csv_output_info->m_last_inserted_tiledb_position + csv_output_info->m_tiledb_anchor_interval;
            //Far off from last inserted position, insert anchor cells
            if(position > curr_anchor_position + csv_output_info->m_tiledb_anchor_interval)
            {
                csv_output_info->m_print_tiledb_anchor_cell = 1;
                for(;curr_anchor_position < position;curr_anchor_position+=csv_output_info->m_tiledb_anchor_interval)
                {
                    csv_output_info->m_curr_tiledb_anchor_position = curr_anchor_position;
                    for(i=0;i<bcf_hdr_nsamples(out_hdr);++i)
                        write_csv_line(mapping_info, csv_output_info, out_hdr, line, i);
                }
                csv_output_info->m_print_tiledb_anchor_cell = 0;
            }
        }
	for(i=0;i<bcf_hdr_nsamples(out_hdr);++i)
	    write_csv_line(mapping_info, csv_output_info, out_hdr, line, i);
    }
}

enum ArgsIdxEnum
{
  ARGS_IDX_SQLITE_FILE=10000,
  ARGS_IDX_COMPUTE_PLMEDIAN,
  ARGS_IDX_SPIT_RANDOM_POSITIONS,
  ARGS_IDX_NUM_RANDOM_POSITIONS,
  ARGS_IDX_NUM_LINES_IN_VCF,
  ARGS_IDX_QUERY_POSITIONS_FILE,
  ARGS_IDX_PROFILE_GVCF_INTERVALS,
  ARGS_IDX_TILEDB_OUTPUT_FORMAT,
  ARGS_IDX_TILEDB_OVERRIDE_SAMPLE_NAME,
  ARGS_IDX_TILEDB_ANCHOR_INTERVAL,
  ARGS_IDX_TILEDB_TREAT_DELETIONS_AS_INTERVALS,
  ARGS_IDX_TAG
};

#define FT_
int main_vcfview(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;
    int targets_is_file = 0, regions_is_file = 0;
    random_sampling_struct* sampling_info = &(args->m_csv_output_info.m_sampling_info);
    gvcf_stat_struct* profile_intervals = &(args->m_csv_output_info.m_profile_intervals);
    static struct option loptions[] =
    {
        {"genotype",1,0,'g'},
        {"compression-level",1,0,'l'},
        {"header-only",0,0,'h'},
        {"no-header",0,0,'H'},
        {"exclude",1,0,'e'},
        {"include",1,0,'i'},
        {"trim-alt-alleles",0,0,'a'},
        {"no-update",0,0,'I'},
        {"drop-genotypes",0,0,'G'},
        {"private",0,0,'x'},
        {"exclude-private",0,0,'X'},
        {"uncalled",0,0,'u'},
        {"exclude-uncalled",0,0,'U'},
        {"apply-filters",1,0,'f'},
        {"known",0,0,'k'},
        {"novel",0,0,'n'},
        {"min-alleles",1,0,'m'},
        {"max-alleles",1,0,'M'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"force-samples",0,0,1},
        {"output-type",1,0,'O'},
        {"output-file",1,0,'o'},
        {"types",1,0,'v'},
        {"exclude-types",1,0,'V'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"min-ac",1,0,'c'},
        {"max-ac",1,0,'C'},
        {"min-af",1,0,'q'},
        {"max-af",1,0,'Q'},
        {"phased",0,0,'p'},
        {"exclude-phased",0,0,'P'},
        {"sqlite",1,0,ARGS_IDX_SQLITE_FILE},
        {"PLmedian",1,0,ARGS_IDX_COMPUTE_PLMEDIAN},
        {"spit-random-positions",1,0,ARGS_IDX_SPIT_RANDOM_POSITIONS},
        {"num-random-positions",1,0,ARGS_IDX_NUM_RANDOM_POSITIONS},
        {"num-lines",1,0,ARGS_IDX_NUM_LINES_IN_VCF},
        {"query-positions-file",1,0,ARGS_IDX_QUERY_POSITIONS_FILE},
	{"profile-gvcf-intervals",0,0,ARGS_IDX_PROFILE_GVCF_INTERVALS},
        {"tiledb-output-format",1,0,ARGS_IDX_TILEDB_OUTPUT_FORMAT},
        {"tiledb-override-sample-name",1,0,ARGS_IDX_TILEDB_OVERRIDE_SAMPLE_NAME},
        {"tiledb-anchor-interval",1,0,ARGS_IDX_TILEDB_ANCHOR_INTERVAL},
        {"tiledb-treat-deletions-as-intervals",0,0,ARGS_IDX_TILEDB_TREAT_DELETIONS_AS_INTERVALS},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "l:t:T:r:R:o:O:s:S:Gf:knv:V:m:M:auUhHc:C:Ii:e:xXpPq:Q:g:",loptions,NULL)) >= 0)
    {
        char allele_type[8] = "nref";
        switch (c)
        {
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    case 't':
                              args->output_type = FT_TILEDB_CSV;
                              args->print_header = 0;
                              break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'l':
                args->clevel = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --compression-level %s\n", optarg);
                args->output_type |= FT_GZ; 
                break;
            case 'o': args->fn_out = optarg; break;
            case 'H': args->print_header = 0; break;
            case 'h': args->header_only = 1; break;

            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;

            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case  1 : args->force_samples = 1; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;

            case 'f': args->m_apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm':
                args->min_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-alleles %s\n", optarg);
                break;
            case 'M': 
                args->max_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --max-alleles %s\n", optarg);
                break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c':
            {
                args->min_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->min_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->min_ac)!=1 )
                    error("Error: Could not parse --min-ac %s\n", optarg);
                set_allele_type(&args->min_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'C':
            {
                args->max_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->max_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->max_ac)!=1 )
                    error("Error: Could not parse --max-ac %s\n", optarg);
                set_allele_type(&args->max_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'q':
            {
                args->min_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->min_af, allele_type)!=2 && sscanf(optarg,"%f",&args->min_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->min_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'Q':
            {
                args->max_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->max_af, allele_type)!=2 && sscanf(optarg,"%f",&args->max_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->max_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }

            case 'x': args->private_vars |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'X': args->private_vars |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'u': args->uncalled |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'U': args->uncalled |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'p': args->phased |= FLT_INCLUDE; break; // phased
            case 'P': args->phased |= FLT_EXCLUDE; break; // exclude-phased
            case 'g':
            {
                if ( !strcasecmp(optarg,"hom") ) args->gt_type = GT_NEED_HOM;
                else if ( !strcasecmp(optarg,"het") ) args->gt_type = GT_NEED_HET;
                else if ( !strcasecmp(optarg,"miss") ) args->gt_type = GT_NEED_MISSING;
                else if ( !strcasecmp(optarg,"^hom") ) args->gt_type = GT_NO_HOM;
                else if ( !strcasecmp(optarg,"^het") ) args->gt_type = GT_NO_HET;
                else if ( !strcasecmp(optarg,"^miss") ) args->gt_type = GT_NO_MISSING;
                else error("The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"%s\".\n", optarg);
                break;
            }
            case ARGS_IDX_SQLITE_FILE:
                strcpy(args->m_mapping_info.sqlite_file, optarg);
                break;
            case ARGS_IDX_COMPUTE_PLMEDIAN:
		args->m_plmedian_info.m_output_fptr = fopen(optarg, "w");
		args->print_header = 0;
		break;
	    case ARGS_IDX_SPIT_RANDOM_POSITIONS:
		sampling_info->m_spit_random_positions = fopen(optarg, "w");
		assert(sampling_info->m_spit_random_positions);
		args->print_header = 0;
		break;
	    case ARGS_IDX_NUM_RANDOM_POSITIONS:
		sampling_info->m_num_random_positions = strtoll(optarg, 0, 10);
		break;
	    case ARGS_IDX_NUM_LINES_IN_VCF:
		sampling_info->m_num_lines_in_vcf = strtoll(optarg, 0, 10);
		break;
	    case ARGS_IDX_QUERY_POSITIONS_FILE:
		args->m_query_positions_file = optarg;
		regions_is_file = 0;
		break;
	    case ARGS_IDX_PROFILE_GVCF_INTERVALS:
		profile_intervals->m_do_profiling = 1;
		break;
            case ARGS_IDX_TILEDB_OUTPUT_FORMAT:
                args->m_mapping_info.m_tiledb_output_version = string_to_tiledb_output_version(optarg);
                break;
            case ARGS_IDX_TILEDB_OVERRIDE_SAMPLE_NAME:
                args->m_mapping_info.m_tiledb_override_sample_name = strdup(optarg);
                break;
            case ARGS_IDX_TILEDB_ANCHOR_INTERVAL:
                args->m_csv_output_info.m_tiledb_anchor_interval = strtoull(optarg, 0, 10);
                break;
            case ARGS_IDX_TILEDB_TREAT_DELETIONS_AS_INTERVALS:
                args->m_csv_output_info.m_treat_deletions_as_intervals = 1;
                break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( args->private_vars > FLT_EXCLUDE ) error("Only one of -x or -X can be given.\n");
    if ( args->uncalled > FLT_EXCLUDE ) error("Only one of -u or -U can be given.\n");
    if ( args->phased > FLT_EXCLUDE ) error("Only one of -p or -P can be given.\n");

    if ( args->sample_names && args->update_info) args->calc_ac = 1;

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    FILE* query_fptr = args->m_query_positions_file ? fopen(args->m_query_positions_file, "r") : 0;
    assert(args->m_query_positions_file == 0 || query_fptr);
    char* query_line = 0;
    size_t query_line_size = 0;
    args->m_first_iteration = 1;
    do
    {
	if(query_fptr)
	{
	    getline(&query_line, &query_line_size, query_fptr);
	    if(feof(query_fptr))
		break;
	    if(query_line[strlen(query_line)-1] == '\n')
		query_line[strlen(query_line)-1] = '\0';
	    args->regions_list = query_line;
	}
	args->files   = bcf_sr_init();
	if(args->m_apply_filters)
	    args->files->apply_filters = args->m_apply_filters;
	// read in the regions from the command line
	if ( args->regions_list )
	{
	    if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
		error("Failed to read the regions: %s\n", args->regions_list);
	}
	else if ( optind+1 < argc )
	{
	    int i;
	    kstring_t tmp = {0,0,0};
	    kputs(argv[optind+1],&tmp);
	    for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
	    if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
		error("Failed to read the regions: %s\n", tmp.s);
	    free(tmp.s);
	}
	if ( args->targets_list )
	{
	    if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
		error("Failed to read the targets: %s\n", args->targets_list);
	}

	if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

	init_data(args);

	bcf_hdr_t *out_hdr = args->hnull ? args->hnull : (args->hsub ? args->hsub : args->hdr);
	if(args->output_type == FT_TILEDB_CSV && args->m_mapping_info.input_sample_idx_2_global_idx == 0)
        {
            initialize_samples_contigs_and_fields_idx(&(args->m_mapping_info), out_hdr);
            initialize_csv_output_info(&(args->m_csv_output_info), args->fn_out ? fopen(args->fn_out, "w") : stdout, 16384, 0);
        }
	if (args->print_header)
	    bcf_hdr_write(args->out, out_hdr);
	else if ( args->output_type & FT_BCF )
	    error("BCF output requires header, cannot proceed with -H\n");
	if(sampling_info->m_spit_random_positions)
	{
	    srand(0x55301b7d);	//random seed
	    assert(sampling_info->m_num_random_positions > 0);
	    assert(sampling_info->m_num_lines_in_vcf > 0);
	    sampling_info->m_sampling_limit = RAND_MAX;
	    sampling_info->m_sampling_limit = 
		(sampling_info->m_sampling_limit * sampling_info->m_num_random_positions)/sampling_info->m_num_lines_in_vcf;
	    assert(sampling_info->m_sampling_limit <= RAND_MAX);
	}
	if (!args->header_only)
	{
	    while ( bcf_sr_next_line(args->files) )
	    {
		bcf1_t *line = args->files->readers[0].buffer[0];
		if ( line->errcode && out_hdr!=args->hdr ) error("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
		if ( subset_vcf(args, line) )
		{
		    if(args->output_type & FT_TILEDB_CSV)	//CSV format
		    {
			if(args->m_plmedian_info.m_output_fptr)	//Query PLmedians
			{
			    int veclen = compute_PLmedian(out_hdr, line,
				    &(args->m_plmedian_info.m_median_result), &(args->m_plmedian_info.m_median_result_len),
				    &(args->m_plmedian_info.m_buffer), &(args->m_plmedian_info.m_buffer_len),
				    &(args->m_plmedian_info.m_reorg_buffer), &(args->m_plmedian_info.m_reorg_buffer_len));
			    print_PLmedian(args->m_plmedian_info.m_output_fptr, out_hdr, line,
				    args->m_plmedian_info.m_median_result, veclen);
			}
			else
			    write_csv_for_one_VCF_line(&(args->m_mapping_info), &(args->m_csv_output_info), out_hdr, line);
		    }
		    else	//VCF/BCF format
		    {
			if(sampling_info->m_spit_random_positions)
			{
			    int random_value = rand();
			    if(random_value < sampling_info->m_sampling_limit)
				fprintf(sampling_info->m_spit_random_positions,"%s\t%d\n",bcf_hdr_id2name(out_hdr,line->rid),line->pos+1); //since line->pos is 0 based, but VCF expects 1 based
			}
			else
			    bcf_write1(args->out, out_hdr, line);
		    }
		}
	    }
	}
	if(!(args->output_type & FT_TILEDB_CSV))
	    hts_close(args->out);
	if(sampling_info->m_spit_random_positions)
	    fclose(sampling_info->m_spit_random_positions);
	bcf_sr_destroy(args->files);
	args->m_first_iteration = 0;
    } while(args->m_query_positions_file && !feof(query_fptr));
    if(profile_intervals->m_do_profiling)
    {
	/*printf("mean valid interval length,mean invalid interval length,,max valid interval length,max invalid interval length,,std-dev valid interval length,std-dev invalid interval length\n");*/
	double mean_valid_interval_length = ((double)profile_intervals->m_num_valid_positions)/profile_intervals->m_num_valid_intervals;
	double stddev_valid_interval_length = sqrt(((double)profile_intervals->m_sum_sq_valid_interval_length)/profile_intervals->m_num_valid_intervals - mean_valid_interval_length*mean_valid_interval_length);
	double mean_invalid_interval_length = ((double)profile_intervals->m_num_invalid_positions)/profile_intervals->m_num_invalid_intervals;
	double stddev_invalid_interval_length = sqrt(((double)profile_intervals->m_sum_sq_invalid_interval_length)/profile_intervals->m_num_invalid_intervals - mean_invalid_interval_length*mean_invalid_interval_length);
	printf("%.2lf,%2.lf,,%"PRIu64",%"PRIu64",,%.2lf,%.2lf\n",mean_valid_interval_length, mean_invalid_interval_length,
		profile_intervals->m_max_valid_interval_length, profile_intervals->m_max_invalid_interval_length,
		stddev_valid_interval_length,stddev_invalid_interval_length);

	/*printf("Max genomic positions: %" PRIu64 "\nNum positions with valid genomic data: %"PRIu64"\nNum valid intervals: %"PRIu64"\n",*/
	/*profile_intervals->m_max_position, profile_intervals->m_num_valid_positions, profile_intervals->m_num_valid_intervals);*/
	/*printf("Num positions with no valid genomic data: %"PRIu64"\nNum invalid intervals: %"PRIu64"\n",*/
	/*profile_intervals->m_num_invalid_positions, profile_intervals->m_num_invalid_intervals);*/
    }
    destroy_data(args);
    destroy_PLmedian(&(args->m_plmedian_info));
    if(query_line)
	free(query_line);
    if(query_fptr)
	fclose(query_fptr);
    free(args);
    return 0;
}

