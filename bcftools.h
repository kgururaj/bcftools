/*  bcftools.h -- utility function declarations.

    Copyright (C) 2013 Genome Research Ltd.

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

#ifndef BCFTOOLS_H
#define BCFTOOLS_H

#include <stdarg.h>
#include <htslib/vcf.h>
#include <math.h>

#define FT_GZ 1
#define FT_VCF 2
#define FT_VCF_GZ (FT_GZ|FT_VCF)
#define FT_BCF (1<<2)
#define FT_BCF_GZ (FT_GZ|FT_BCF)
#define FT_STDIN (1<<3)
#define FT_TILEDB_CSV (1<<4)

#ifdef __cplusplus
extern "C"
{
#endif

  char *bcftools_version(void);
  void error(const char *format, ...);
  void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd);
  const char *hts_bcf_wmode(int file_type);

  void *smalloc(size_t size);     // safe malloc

  static inline char gt2iupac(char a, char b)
  {
    static const char iupac[4][4] = { {'A','M','R','W'},{'M','C','S','Y'},{'R','S','G','K'},{'W','Y','K','T'} };
    if ( a>='a' ) a -= 'a' - 'A';
    if ( b>='a' ) b -= 'a' - 'A';
    if ( a=='A' ) a = 0;
    else if ( a=='C' ) a = 1;
    else if ( a=='G' ) a = 2;
    else if ( a=='T' ) a = 3;
    else return 'N';
    if ( b=='A' ) b = 0;
    else if ( b=='C' ) b = 1;
    else if ( b=='G' ) b = 2;
    else if ( b=='T' ) b = 3;
    else return 'N';
    return iupac[(int)a][(int)b];
  }

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

  //Median of PL functions
  int compute_PLmedian(bcf_hdr_t* hdr, bcf1_t* line, int** median_result, int* median_result_len,
      int** buffer, int* buffer_len,
      int** reorg_buffer, int* reorg_buffer_len);
  void print_PLmedian(FILE* fptr, bcf_hdr_t* hdr, bcf1_t* line, int* median_result, int median_result_len);
  void destroy_PLmedian(plmedian_struct* info);

  //TileDB functions and structures
#include <sqlite3.h>
  typedef struct
  {
    unsigned m_offset;
    unsigned m_size;
    char* m_buffer;
    int64_t m_global_sample_idx;
    int64_t m_global_start_position;
    int64_t m_global_end_position;
  } buffer_wrapper;

#define ALLOCATE_SPACE_IN_BUFFER(buffer, bytes) \
  buffer->m_offset += bytes; \
  if(buffer->m_offset >= buffer->m_size) \
      error("TileDB CSV output buffer has overflown :(\n");

  enum TileDBOutputVersion
  {
    TILEDB_OUTPUT_CSV_V0=0,
    TILEDB_OUTPUT_CSV_V1,
    TILEDB_OUTPUT_BINARY_V1,
  };
  typedef void (*TileDBPrinterTy) (buffer_wrapper* buffer, unsigned bcf_ht_type, const char* fmt, ...);
  typedef struct
  {
    unsigned m_tiledb_output_version;
    //Pointer to printer function
    TileDBPrinterTy m_tiledb_printer;
    //List of field names
    char** m_field_names;
    char** m_contig_names;
    uint64_t* m_contig_lengths;
    //Mappings
    //Value in ith position specifies global sample idx in SQLite DB for ith
    //sample in the query array
    int64_t* input_sample_idx_2_global_idx;
    int64_t* input_field_idx_2_global_idx;
    int64_t* input_contig_idx_2_global_idx;
    int64_t* input_contig_idx_2_offset;
    char sqlite_file[1024];
    sqlite3* db;
  }sqlite_mappings_struct;

    typedef struct
    {
      int m_input_sample_idx;
      int m_reader_idx;
      int64_t m_global_sample_idx;
    }global_samples_struct;

  typedef struct
  {
    int m_max_contig_idx;
    int m_max_sample_idx;
    int m_current_contig_idx;
    sqlite_mappings_struct m_merged_sqlite_mapping_info;
    int64_t m_last_global_position;
    char* m_output_directory;
    FILE** m_contig_output_csv; //1 file per contig
    char* m_sqlite_file;
    sqlite3* m_sqlite_db;
    global_samples_struct* m_globally_sorted_sample_info;
    global_samples_struct** m_current_position_samples_array;
    int64_t m_num_current_position_samples;
  }tiledb_struct;

  typedef struct
  {
    FILE* m_spit_random_positions;	//0 in calloc
    int64_t m_num_random_positions;
    int64_t m_num_lines_in_vcf;
    int64_t m_sampling_limit;
  }random_sampling_struct;

#define SHOULD_DO_SAMPLING(X) ((X).m_spit_random_positions)

  typedef struct
  {
    uint8_t m_do_profiling;
    uint64_t m_num_valid_positions;
    uint64_t m_sum_sq_valid_interval_length;
    uint64_t m_max_valid_interval_length;
    uint64_t m_num_valid_intervals;
    uint64_t m_max_position;
    uint64_t m_num_invalid_positions;
    uint64_t m_sum_sq_invalid_interval_length;
    uint64_t m_max_invalid_interval_length;
    uint64_t m_num_invalid_intervals;
    uint64_t m_last_valid_interval_end;
  }gvcf_stat_struct;

#define SHOULD_DO_PROFILING(X) ((X).m_do_profiling)

  //Big bag of output control items
  typedef struct
  {
    uint8_t m_skip_coordinates;
    buffer_wrapper m_csv_out_buffer;
    FILE* m_csv_out_fptr; //non-0 if output to file/stream
    gvcf_stat_struct m_profile_intervals; //non-0 if profiling info to be collected
    random_sampling_struct m_sampling_info;//non-0 if sampling needs to be done
    char* m_htslib_buffer;
    unsigned m_htslib_buffer_size; 
  } csv_output_struct;

  void open_sqlite3_db(const char* sqlite_file, sqlite3** db);
  void initialize_samples_contigs_and_fields_idx(sqlite_mappings_struct* mapping_info, const bcf_hdr_t* hdr);
  /*
   * Allocate sqlite_mappings_struct and open SQLite file
   * @return Pointer to sqlite_mappings_struct with opened sqlite_file
   * */
  void* allocate_sqlite3_mapping(const char* sqlite_file);
  /*
   * Query sqlite DB function
   * @param info_ptr: sqlite_mappings_struct* returned by allocate_sqlite3_mapping() function. Members of this structure
   * are modified by this function
   * @param n_samples: #samples being queried
   * @param sample_names: array of sample names char**
   * @param contig_lengths: if the query is expected to add new contigs to the
   * SQLite DB, then this arg should point to an array of contig lengths. Can
   * be NULL, if no updates are expected to be performed.
   * @return Pointer to array containing result of query
   * */
  int64_t* query_samples_idx(void* info_ptr, int n_samples, const char* const* sample_names);
  int64_t* query_contigs_offset(void* info_ptr, int n_contigs, const char* const* contig_names, const uint64_t* contig_lengths);
  int64_t* query_fields_idx(void* info_ptr, int n_fields, const char* const* field_names);
  /*
   * Free members of sqlite_mappings_struct info_ptr
   */
  void free_sqlite3_data(void* info_ptr);
  void initialize_csv_output_info(csv_output_struct* ptr, FILE* output_fptr, unsigned buffer_size, uint8_t reset);
  void free_csv_output_info(csv_output_struct* info);
  //Write CSV line for a single sample (input_sample_idx)
  int write_csv_line(sqlite_mappings_struct* mapping_info, csv_output_struct* csv_output_info,
      bcf_hdr_t* out_hdr, bcf1_t* line, int input_sample_idx);
  //Write CSV lines for all samples in a given line in the VCF
  void write_csv_for_one_VCF_line(sqlite_mappings_struct* mapping_info, csv_output_struct* csv_output_info,
      bcf_hdr_t* out_hdr, bcf1_t* line);

  static inline double phred_score(double prob)
  {
    if ( prob==0 ) return 99;
    prob = -4.3429*log(prob);
    return prob>99 ? 99 : prob;
  }

  unsigned string_to_tiledb_output_version(const char* string);
#ifdef __cplusplus
}
#endif

#endif
