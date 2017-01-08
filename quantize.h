#pragma once

#define RGB_RED       0
#define RGB_GREEN     1
#define RGB_BLUE      2

#define MAXJSAMPLE        255
#define CENTERJSAMPLE     128
#define BITS_IN_JSAMPLE   8
#define GETJSAMPLE(value) ((int) (value))

#define RIGHT_SHIFT(x,shft) ((x) >> (shft))

#define MAXNUMCOLORS  (MAXJSAMPLE+1) /* maximum size of colormap */

/* These will do the right thing for either R,G,B or B,G,R color order,
* but you may not like the results for other color orders.
*/
#define HIST_C0_BITS  5     /* bits of precision in R/B histogram */
#define HIST_C1_BITS  6     /* bits of precision in G histogram */
#define HIST_C2_BITS  5     /* bits of precision in B/R histogram */

/* Number of elements along histogram axes. */
#define HIST_C0_ELEMS  (1<<HIST_C0_BITS)
#define HIST_C1_ELEMS  (1<<HIST_C1_BITS)
#define HIST_C2_ELEMS  (1<<HIST_C2_BITS)

/* These are the amounts to shift an input value to get a histogram index. */
#define C0_SHIFT  (BITS_IN_JSAMPLE-HIST_C0_BITS)
#define C1_SHIFT  (BITS_IN_JSAMPLE-HIST_C1_BITS)
#define C2_SHIFT  (BITS_IN_JSAMPLE-HIST_C2_BITS)

/* log2(histogram cells in update box) for each axis; this can be adjusted */
#define BOX_C0_LOG  (HIST_C0_BITS-3)
#define BOX_C1_LOG  (HIST_C1_BITS-3)
#define BOX_C2_LOG  (HIST_C2_BITS-3)

#define BOX_C0_ELEMS  (1<<BOX_C0_LOG) /* # of hist cells in update box */
#define BOX_C1_ELEMS  (1<<BOX_C1_LOG)
#define BOX_C2_ELEMS  (1<<BOX_C2_LOG)

#define BOX_C0_SHIFT  (C0_SHIFT + BOX_C0_LOG)
#define BOX_C1_SHIFT  (C1_SHIFT + BOX_C1_LOG)
#define BOX_C2_SHIFT  (C2_SHIFT + BOX_C2_LOG)

#define R_SCALE 2       /* scale R distances by this much */
#define G_SCALE 3       /* scale G distances by this much */
#define B_SCALE 1       /* and B by this much */

#if RGB_RED == 0
#define C0_SCALE R_SCALE
#endif
#if RGB_BLUE == 0
#define C0_SCALE B_SCALE
#endif
#if RGB_GREEN == 1
#define C1_SCALE G_SCALE
#endif
#if RGB_RED == 2
#define C2_SCALE R_SCALE
#endif
#if RGB_BLUE == 2
#define C2_SCALE B_SCALE
#endif

typedef JSAMPROW *JSAMPARRAY;
typedef unsigned int JDIMENSION;

typedef wxUint16 histcell;    /* histogram cell; prefer an unsigned type */

typedef histcell  * histptr;    /* for pointers to histogram cells */

typedef histcell hist1d[HIST_C2_ELEMS]; /* typedefs for the array */
typedef hist1d  * hist2d;   /* type for the 2nd-level pointers */
typedef hist2d * hist3d;    /* type for top-level pointer */

typedef struct {
	void *cquantize;
	JDIMENSION output_width;
	JSAMPARRAY colormap;
	int actual_number_of_colors;
	int desired_number_of_colors;
	JSAMPLE *sample_range_limit, *srl_orig;
} j_decompress;

#if defined(__WINDOWS__)
#define JMETHOD(type,methodname,arglist)  type (__cdecl methodname) arglist
#else
#define JMETHOD(type,methodname,arglist)  type (methodname) arglist
#endif

typedef j_decompress *j_decompress_ptr;
struct jpeg_color_quantizer {
	JMETHOD(void, start_pass, (j_decompress_ptr cinfo, bool is_pre_scan));
	JMETHOD(void, color_quantize, (j_decompress_ptr cinfo,
		JSAMPARRAY input_buf, JSAMPARRAY output_buf,
		int num_rows));
	JMETHOD(void, finish_pass, (j_decompress_ptr cinfo));
	JMETHOD(void, new_color_map, (j_decompress_ptr cinfo));
};


typedef wxInt16 FSERROR;      /* 16 bits should be enough */
typedef int LOCFSERROR;     /* use 'int' for calculation temps */

typedef FSERROR  *FSERRPTR; /* pointer to error array (in  storage!) */

typedef struct {

	struct {
		void(*finish_pass)(j_decompress_ptr);
		void(*color_quantize)(j_decompress_ptr, JSAMPARRAY, JSAMPARRAY, int);
		void(*start_pass)(j_decompress_ptr, bool);
		void(*new_color_map)(j_decompress_ptr);
	} pub;

	/* Space for the eventually created colormap is stashed here */
	JSAMPARRAY sv_colormap;   /* colormap allocated at init time */
	int desired;          /* desired # of colors = size of colormap */

						  /* Variables for accumulating image statistics */
	hist3d histogram;     /* pointer to the histogram */

	bool needs_zeroed;        /* true if next pass must zero histogram */

							  /* Variables for Floyd-Steinberg dithering */
	FSERRPTR fserrors;        /* accumulated errors */
	bool on_odd_row;      /* flag to remember which row we are on */
	int * error_limiter;      /* table for clamping the applied error */
} my_cquantizer;

typedef my_cquantizer * my_cquantize_ptr;

void
fill_inverse_cmap(j_decompress_ptr cinfo, int c0, int c1, int c2);
void
pass2_fs_dither(j_decompress_ptr cinfo,
	JSAMPARRAY input_buf, JSAMPARRAY output_buf, int num_rows);
