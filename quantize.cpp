#include "stdafx.h"

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

static int
find_nearby_colors(j_decompress_ptr cinfo, int minc0, int minc1, int minc2,
	JSAMPLE colorlist[])
	/* Locate the colormap entries close enough to an update box to be candidates
	* for the nearest entry to some cell(s) in the update box.  The update box
	* is specified by the center coordinates of its first cell.  The number of
	* candidate colormap entries is returned, and their colormap indexes are
	* placed in colorlist[].
	* This routine uses Heckbert's "locally sorted search" criterion to select
	* the colors that need further consideration.
	*/
{
	int numcolors = cinfo->actual_number_of_colors;
	int maxc0, maxc1, maxc2;
	int centerc0, centerc1, centerc2;
	int i, x, ncolors;
	wxInt32 minmaxdist, min_dist, max_dist, tdist;
	wxInt32 mindist[MAXNUMCOLORS];  /* min distance to colormap entry i */

									/* Compute true coordinates of update box's upper corner and center.
									* Actually we compute the coordinates of the center of the upper-corner
									* histogram cell, which are the upper bounds of the volume we care about.
									* Note that since ">>" rounds down, the "center" values may be closer to
									* min than to max; hence comparisons to them must be "<=", not "<".
									*/
	maxc0 = minc0 + ((1 << BOX_C0_SHIFT) - (1 << C0_SHIFT));
	centerc0 = (minc0 + maxc0) >> 1;
	maxc1 = minc1 + ((1 << BOX_C1_SHIFT) - (1 << C1_SHIFT));
	centerc1 = (minc1 + maxc1) >> 1;
	maxc2 = minc2 + ((1 << BOX_C2_SHIFT) - (1 << C2_SHIFT));
	centerc2 = (minc2 + maxc2) >> 1;

	/* For each color in colormap, find:
	*  1. its minimum squared-distance to any point in the update box
	*     (zero if color is within update box);
	*  2. its maximum squared-distance to any point in the update box.
	* Both of these can be found by considering only the corners of the box.
	* We save the minimum distance for each color in mindist[];
	* only the smallest maximum distance is of interest.
	*/
	minmaxdist = 0x7FFFFFFFL;

	for (i = 0; i < numcolors; i++) {
		/* We compute the squared-c0-distance term, then add in the other two. */
		x = GETJSAMPLE(cinfo->colormap[0][i]);
		if (x < minc0) {
			tdist = (x - minc0) * C0_SCALE;
			min_dist = tdist*tdist;
			tdist = (x - maxc0) * C0_SCALE;
			max_dist = tdist*tdist;
		}
		else if (x > maxc0) {
			tdist = (x - maxc0) * C0_SCALE;
			min_dist = tdist*tdist;
			tdist = (x - minc0) * C0_SCALE;
			max_dist = tdist*tdist;
		}
		else {
			/* within cell range so no contribution to min_dist */
			min_dist = 0;
			if (x <= centerc0) {
				tdist = (x - maxc0) * C0_SCALE;
				max_dist = tdist*tdist;
			}
			else {
				tdist = (x - minc0) * C0_SCALE;
				max_dist = tdist*tdist;
			}
		}

		x = GETJSAMPLE(cinfo->colormap[1][i]);
		if (x < minc1) {
			tdist = (x - minc1) * C1_SCALE;
			min_dist += tdist*tdist;
			tdist = (x - maxc1) * C1_SCALE;
			max_dist += tdist*tdist;
		}
		else if (x > maxc1) {
			tdist = (x - maxc1) * C1_SCALE;
			min_dist += tdist*tdist;
			tdist = (x - minc1) * C1_SCALE;
			max_dist += tdist*tdist;
		}
		else {
			/* within cell range so no contribution to min_dist */
			if (x <= centerc1) {
				tdist = (x - maxc1) * C1_SCALE;
				max_dist += tdist*tdist;
			}
			else {
				tdist = (x - minc1) * C1_SCALE;
				max_dist += tdist*tdist;
			}
		}

		x = GETJSAMPLE(cinfo->colormap[2][i]);
		if (x < minc2) {
			tdist = (x - minc2) * C2_SCALE;
			min_dist += tdist*tdist;
			tdist = (x - maxc2) * C2_SCALE;
			max_dist += tdist*tdist;
		}
		else if (x > maxc2) {
			tdist = (x - maxc2) * C2_SCALE;
			min_dist += tdist*tdist;
			tdist = (x - minc2) * C2_SCALE;
			max_dist += tdist*tdist;
		}
		else {
			/* within cell range so no contribution to min_dist */
			if (x <= centerc2) {
				tdist = (x - maxc2) * C2_SCALE;
				max_dist += tdist*tdist;
			}
			else {
				tdist = (x - minc2) * C2_SCALE;
				max_dist += tdist*tdist;
			}
		}

		mindist[i] = min_dist;  /* save away the results */
		if (max_dist < minmaxdist)
			minmaxdist = max_dist;
	}

	/* Now we know that no cell in the update box is more than minmaxdist
	* away from some colormap entry.  Therefore, only colors that are
	* within minmaxdist of some part of the box need be considered.
	*/
	ncolors = 0;
	for (i = 0; i < numcolors; i++) {
		if (mindist[i] <= minmaxdist)
			colorlist[ncolors++] = (JSAMPLE)i;
	}
	return ncolors;
}


static void
find_best_colors(j_decompress_ptr cinfo, int minc0, int minc1, int minc2,
	int numcolors, JSAMPLE colorlist[], JSAMPLE bestcolor[])
	/* Find the closest colormap entry for each cell in the update box,
	* given the list of candidate colors prepared by find_nearby_colors.
	* Return the indexes of the closest entries in the bestcolor[] array.
	* This routine uses Thomas' incremental distance calculation method to
	* find the distance from a colormap entry to successive cells in the box.
	*/
{
	int ic0, ic1, ic2;
	int i, icolor;
	wxInt32 * bptr;    /* pointer into bestdist[] array */
	JSAMPLE * cptr;       /* pointer into bestcolor[] array */
	wxInt32 dist0, dist1;       /* initial distance values */
	wxInt32 dist2;     /* current distance in inner loop */
	wxInt32 xx0, xx1;       /* distance increments */
	wxInt32 xx2;
	wxInt32 inc0, inc1, inc2;   /* initial values for increments */
								/* This array holds the distance to the nearest-so-far color for each cell */
	wxInt32 bestdist[BOX_C0_ELEMS * BOX_C1_ELEMS * BOX_C2_ELEMS];

	/* Initialize best-distance for each cell of the update box */
	bptr = bestdist;
	for (i = BOX_C0_ELEMS*BOX_C1_ELEMS*BOX_C2_ELEMS - 1; i >= 0; i--)
		*bptr++ = 0x7FFFFFFFL;

	/* For each color selected by find_nearby_colors,
	* compute its distance to the center of each cell in the box.
	* If that's less than best-so-far, update best distance and color number.
	*/

	/* Nominal steps between cell centers ("x" in Thomas article) */
#define STEP_C0  ((1 << C0_SHIFT) * C0_SCALE)
#define STEP_C1  ((1 << C1_SHIFT) * C1_SCALE)
#define STEP_C2  ((1 << C2_SHIFT) * C2_SCALE)

	for (i = 0; i < numcolors; i++) {
		icolor = GETJSAMPLE(colorlist[i]);
		/* Compute (square of) distance from minc0/c1/c2 to this color */
		inc0 = (minc0 - GETJSAMPLE(cinfo->colormap[0][icolor])) * C0_SCALE;
		dist0 = inc0*inc0;
		inc1 = (minc1 - GETJSAMPLE(cinfo->colormap[1][icolor])) * C1_SCALE;
		dist0 += inc1*inc1;
		inc2 = (minc2 - GETJSAMPLE(cinfo->colormap[2][icolor])) * C2_SCALE;
		dist0 += inc2*inc2;
		/* Form the initial difference increments */
		inc0 = inc0 * (2 * STEP_C0) + STEP_C0 * STEP_C0;
		inc1 = inc1 * (2 * STEP_C1) + STEP_C1 * STEP_C1;
		inc2 = inc2 * (2 * STEP_C2) + STEP_C2 * STEP_C2;
		/* Now loop over all cells in box, updating distance per Thomas method */
		bptr = bestdist;
		cptr = bestcolor;
		xx0 = inc0;
		for (ic0 = BOX_C0_ELEMS - 1; ic0 >= 0; ic0--) {
			dist1 = dist0;
			xx1 = inc1;
			for (ic1 = BOX_C1_ELEMS - 1; ic1 >= 0; ic1--) {
				dist2 = dist1;
				xx2 = inc2;
				for (ic2 = BOX_C2_ELEMS - 1; ic2 >= 0; ic2--) {
					if (dist2 < *bptr) {
						*bptr = dist2;
						*cptr = (JSAMPLE)icolor;
					}
					dist2 += xx2;
					xx2 += 2 * STEP_C2 * STEP_C2;
					bptr++;
					cptr++;
				}
				dist1 += xx1;
				xx1 += 2 * STEP_C1 * STEP_C1;
			}
			dist0 += xx0;
			xx0 += 2 * STEP_C0 * STEP_C0;
		}
	}
}

static void
fill_inverse_cmap(j_decompress_ptr cinfo, int c0, int c1, int c2)
/* Fill the inverse-colormap entries in the update box that contains */
/* histogram cell c0/c1/c2.  (Only that one cell MUST be filled, but */
/* we can fill as many others as we wish.) */
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	hist3d histogram = cquantize->histogram;
	int minc0, minc1, minc2;  /* lower left corner of update box */
	int ic0, ic1, ic2;
	JSAMPLE * cptr;  /* pointer into bestcolor[] array */
	histptr cachep;  /* pointer into main cache array */
					 /* This array lists the candidate colormap indexes. */
	JSAMPLE colorlist[MAXNUMCOLORS];
	int numcolors;        /* number of candidate colors */
						  /* This array holds the actually closest colormap index for each cell. */
	JSAMPLE bestcolor[BOX_C0_ELEMS * BOX_C1_ELEMS * BOX_C2_ELEMS];

	/* Convert cell coordinates to update box ID */
	c0 >>= BOX_C0_LOG;
	c1 >>= BOX_C1_LOG;
	c2 >>= BOX_C2_LOG;

	/* Compute true coordinates of update box's origin corner.
	* Actually we compute the coordinates of the center of the corner
	* histogram cell, which are the lower bounds of the volume we care about.
	*/
	minc0 = (c0 << BOX_C0_SHIFT) + ((1 << C0_SHIFT) >> 1);
	minc1 = (c1 << BOX_C1_SHIFT) + ((1 << C1_SHIFT) >> 1);
	minc2 = (c2 << BOX_C2_SHIFT) + ((1 << C2_SHIFT) >> 1);

	/* Determine which colormap entries are close enough to be candidates
	* for the nearest entry to some cell in the update box.
	*/
	numcolors = find_nearby_colors(cinfo, minc0, minc1, minc2, colorlist);

	/* Determine the actually nearest colors. */
	find_best_colors(cinfo, minc0, minc1, minc2, numcolors, colorlist,
		bestcolor);

	/* Save the best color numbers (plus 1) in the main cache array */
	c0 <<= BOX_C0_LOG;        /* convert ID back to base cell indexes */
	c1 <<= BOX_C1_LOG;
	c2 <<= BOX_C2_LOG;
	cptr = bestcolor;
	for (ic0 = 0; ic0 < BOX_C0_ELEMS; ic0++) {
		for (ic1 = 0; ic1 < BOX_C1_ELEMS; ic1++) {
			cachep = &histogram[c0 + ic0][c1 + ic1][c2];
			for (ic2 = 0; ic2 < BOX_C2_ELEMS; ic2++) {
				*cachep++ = (histcell)(GETJSAMPLE(*cptr++) + 1);
			}
		}
	}
}

void
pass2_fs_dither(j_decompress_ptr cinfo,
	JSAMPARRAY input_buf, JSAMPARRAY output_buf, int num_rows)
	/* This version performs Floyd-Steinberg dithering */
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	hist3d histogram = cquantize->histogram;

	// Making the curN variables "volatile" lets the optimizing compiler
	// generate code that generates the same output as when the optimizer
	// is disabled.
#if 0
	volatile
#endif
	LOCFSERROR cur0, cur1, cur2; /* current error or pixel value */
	LOCFSERROR belowerr0, belowerr1, belowerr2; /* error for pixel below cur */
	LOCFSERROR bpreverr0, bpreverr1, bpreverr2; /* error for below/prev col */
	FSERRPTR errorptr;   /* => fserrors[] at column before current */
	JSAMPROW inptr;       /* => current input pixel */
	JSAMPROW outptr;      /* => current output pixel */
	histptr cachep;
	int dir;          /* +1 or -1 depending on direction */
	int dir3;         /* 3*dir, for advancing inptr & errorptr */
	int row;
	JDIMENSION col;
	JDIMENSION width = cinfo->output_width;
	JSAMPLE *range_limit = cinfo->sample_range_limit;
	int *error_limit = cquantize->error_limiter;
	JSAMPROW colormap0 = cinfo->colormap[0];
	JSAMPROW colormap1 = cinfo->colormap[1];
	JSAMPROW colormap2 = cinfo->colormap[2];


	for (row = 0; row < num_rows; row++) {
		inptr = input_buf[row];
		outptr = output_buf[row];
		if (cquantize->on_odd_row) {
			/* work right to left in this row */
			inptr += (width - 1) * 3;   /* so point to rightmost pixel */
			outptr += width - 1;
			dir = -1;
			dir3 = -3;
			errorptr = cquantize->fserrors + (width + 1) * 3; /* => entry after last column */
			cquantize->on_odd_row = false; /* flip for next time */
		}
		else {
			/* work left to right in this row */
			dir = 1;
			dir3 = 3;
			errorptr = cquantize->fserrors; /* => entry before first real column */
			cquantize->on_odd_row = true; /* flip for next time */
		}
		/* Preset error values: no error propagated to first pixel from left */
		cur0 = cur1 = cur2 = 0;
		/* and no error propagated to row below yet */
		belowerr0 = belowerr1 = belowerr2 = 0;
		bpreverr0 = bpreverr1 = bpreverr2 = 0;

		for (col = width; col > 0; col--) {
			/* curN holds the error propagated from the previous pixel on the
			* current line.  Add the error propagated from the previous line
			* to form the complete error correction term for this pixel, and
			* round the error term (which is expressed * 16) to an integer.
			* RIGHT_SHIFT rounds towards minus infinity, so adding 8 is correct
			* for either sign of the error value.
			* Note: errorptr points to *previous* column's array entry.
			*/
			cur0 = RIGHT_SHIFT(cur0 + errorptr[dir3 + 0] + 8, 4);
			cur1 = RIGHT_SHIFT(cur1 + errorptr[dir3 + 1] + 8, 4);
			cur2 = RIGHT_SHIFT(cur2 + errorptr[dir3 + 2] + 8, 4);
			/* Limit the error using transfer function set by init_error_limit.
			* See comments with init_error_limit for rationale.
			*/
			cur0 = error_limit[cur0];
			cur1 = error_limit[cur1];

			if (cur2 > MAXJSAMPLE || cur2 < -MAXJSAMPLE)
				::OutputDebugStringW(L"Hello\n");

			cur2 = error_limit[cur2];
			/* Form pixel value + error, and range-limit to 0..MAXJSAMPLE.
			* The maximum error is +- MAXJSAMPLE (or less with error limiting);
			* this sets the required size of the range_limit array.
			*/
			cur0 += GETJSAMPLE(inptr[0]);
			cur1 += GETJSAMPLE(inptr[1]);
			cur2 += GETJSAMPLE(inptr[2]);
			cur0 = GETJSAMPLE(range_limit[cur0]);
			cur1 = GETJSAMPLE(range_limit[cur1]);
			cur2 = GETJSAMPLE(range_limit[cur2]);
			/* Index into the cache with adjusted pixel value */
			cachep = &histogram[cur0 >> C0_SHIFT][cur1 >> C1_SHIFT][cur2 >> C2_SHIFT];
			/* If we have not seen this color before, find nearest colormap */
			/* entry and update the cache */
			if (*cachep == 0)
				fill_inverse_cmap(cinfo, cur0 >> C0_SHIFT, cur1 >> C1_SHIFT, cur2 >> C2_SHIFT);
			/* Now emit the colormap index for this cell */
			{ int pixcode = *cachep - 1;
			*outptr = (JSAMPLE)pixcode;
			/* Compute representation error for this pixel */
			cur0 -= GETJSAMPLE(colormap0[pixcode]);
			cur1 -= GETJSAMPLE(colormap1[pixcode]);
			cur2 -= GETJSAMPLE(colormap2[pixcode]);
			}
			/* Compute error fractions to be propagated to adjacent pixels.
			* Add these into the running sums, and simultaneously shift the
			* next-line error sums left by 1 column.
			*/
			{ LOCFSERROR bnexterr, delta;

			bnexterr = cur0;    /* Process component 0 */
			delta = cur0 * 2;
			cur0 += delta;      /* form error * 3 */
			errorptr[0] = (FSERROR)(bpreverr0 + cur0);
			cur0 += delta;      /* form error * 5 */
			bpreverr0 = belowerr0 + cur0;
			belowerr0 = bnexterr;
			cur0 += delta;      /* form error * 7 */
			bnexterr = cur1;    /* Process component 1 */
			delta = cur1 * 2;
			cur1 += delta;      /* form error * 3 */
			errorptr[1] = (FSERROR)(bpreverr1 + cur1);
			cur1 += delta;      /* form error * 5 */
			bpreverr1 = belowerr1 + cur1;
			belowerr1 = bnexterr;
			cur1 += delta;      /* form error * 7 */
			bnexterr = cur2;    /* Process component 2 */
			delta = cur2 * 2;
			cur2 += delta;      /* form error * 3 */
			errorptr[2] = (FSERROR)(bpreverr2 + cur2);
			cur2 += delta;      /* form error * 5 */
			bpreverr2 = belowerr2 + cur2;
			belowerr2 = bnexterr;
			cur2 += delta;      /* form error * 7 */
			}
			/* At this point curN contains the 7/16 error value to be propagated
			* to the next pixel on the current line, and all the errors for the
			* next line have been shifted over.  We are therefore ready to move on.
			*/
			inptr += dir3;        /* Advance pixel pointers to next column */
			outptr += dir;
			errorptr += dir3;     /* advance errorptr to current column */
		}
		/* Post-loop cleanup: we must unload the final error values into the
		* final fserrors[] entry.  Note we need not unload belowerrN because
		* it is for the dummy column before or after the actual array.
		*/
		errorptr[0] = (FSERROR)bpreverr0; /* unload prev errs into array */
		errorptr[1] = (FSERROR)bpreverr1;
		errorptr[2] = (FSERROR)bpreverr2;
	}
}

static void
init_error_limit(j_decompress_ptr cinfo)
/* Allocate and fill in the error_limiter table */
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	int * table;
	int in, out;

	table = (int *)malloc((MAXJSAMPLE * 2 + 1) * sizeof(int));
	table += MAXJSAMPLE;      /* so can index -MAXJSAMPLE .. +MAXJSAMPLE */
	cquantize->error_limiter = table;

#define STEPSIZE ((MAXJSAMPLE+1)/16)
	/* Map errors 1:1 up to +- MAXJSAMPLE/16 */
	out = 0;
	for (in = 0; in < STEPSIZE; in++, out++) {
		table[in] = out; table[-in] = -out;
	}
	/* Map errors 1:2 up to +- 3*MAXJSAMPLE/16 */
	for (; in < STEPSIZE * 3; in++, out += (in & 1) ? 0 : 1) {
		table[in] = out; table[-in] = -out;
	}
	/* Clamp the rest to final out value (which is (MAXJSAMPLE+1)/8) */
	for (; in <= MAXJSAMPLE; in++) {
		table[in] = out; table[-in] = -out;
	}
#undef STEPSIZE
}


typedef struct {
	/* The bounds of the box (inclusive); expressed as histogram indexes */
	int c0min, c0max;
	int c1min, c1max;
	int c2min, c2max;
	/* The volume (actually 2-norm) of the box */
	wxInt32 volume;
	/* The number of nonzero histogram cells within this box */
	long colorcount;
} box;

typedef box * boxptr;


boxptr
find_biggest_color_pop(boxptr boxlist, int numboxes)
/* Find the splittable box with the largest color population */
/* Returns NULL if no splittable boxes remain */
{
	boxptr boxp;
	int i;
	long maxc = 0;
	boxptr which = NULL;

	for (i = 0, boxp = boxlist; i < numboxes; i++, boxp++) {
		if (boxp->colorcount > maxc && boxp->volume > 0) {
			which = boxp;
			maxc = boxp->colorcount;
		}
	}
	return which;
}


boxptr
find_biggest_volume(boxptr boxlist, int numboxes)
/* Find the splittable box with the largest (scaled) volume */
/* Returns NULL if no splittable boxes remain */
{
	boxptr boxp;
	int i;
	wxInt32 maxv = 0;
	boxptr which = NULL;

	for (i = 0, boxp = boxlist; i < numboxes; i++, boxp++) {
		if (boxp->volume > maxv) {
			which = boxp;
			maxv = boxp->volume;
		}
	}
	return which;
}

void
update_box(j_decompress_ptr cinfo, boxptr boxp)
/* Shrink the min/max bounds of a box to enclose only nonzero elements, */
/* and recompute its volume and population */
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	hist3d histogram = cquantize->histogram;
	histptr histp;
	int c0, c1, c2;
	int c0min, c0max, c1min, c1max, c2min, c2max;
	wxInt32 dist0, dist1, dist2;
	long ccount;

	c0min = boxp->c0min;  c0max = boxp->c0max;
	c1min = boxp->c1min;  c1max = boxp->c1max;
	c2min = boxp->c2min;  c2max = boxp->c2max;

	if (c0max > c0min)
		for (c0 = c0min; c0 <= c0max; c0++)
			for (c1 = c1min; c1 <= c1max; c1++) {
				histp = &histogram[c0][c1][c2min];
				for (c2 = c2min; c2 <= c2max; c2++)
					if (*histp++ != 0) {
						boxp->c0min = c0min = c0;
						goto have_c0min;
					}
			}
have_c0min:
	if (c0max > c0min)
		for (c0 = c0max; c0 >= c0min; c0--)
			for (c1 = c1min; c1 <= c1max; c1++) {
				histp = &histogram[c0][c1][c2min];
				for (c2 = c2min; c2 <= c2max; c2++)
					if (*histp++ != 0) {
						boxp->c0max = c0max = c0;
						goto have_c0max;
					}
			}
have_c0max:
	if (c1max > c1min)
		for (c1 = c1min; c1 <= c1max; c1++)
			for (c0 = c0min; c0 <= c0max; c0++) {
				histp = &histogram[c0][c1][c2min];
				for (c2 = c2min; c2 <= c2max; c2++)
					if (*histp++ != 0) {
						boxp->c1min = c1min = c1;
						goto have_c1min;
					}
			}
have_c1min:
	if (c1max > c1min)
		for (c1 = c1max; c1 >= c1min; c1--)
			for (c0 = c0min; c0 <= c0max; c0++) {
				histp = &histogram[c0][c1][c2min];
				for (c2 = c2min; c2 <= c2max; c2++)
					if (*histp++ != 0) {
						boxp->c1max = c1max = c1;
						goto have_c1max;
					}
			}
have_c1max:
	if (c2max > c2min)
		for (c2 = c2min; c2 <= c2max; c2++)
			for (c0 = c0min; c0 <= c0max; c0++) {
				histp = &histogram[c0][c1min][c2];
				for (c1 = c1min; c1 <= c1max; c1++, histp += HIST_C2_ELEMS)
					if (*histp != 0) {
						boxp->c2min = c2min = c2;
						goto have_c2min;
					}
			}
have_c2min:
	if (c2max > c2min)
		for (c2 = c2max; c2 >= c2min; c2--)
			for (c0 = c0min; c0 <= c0max; c0++) {
				histp = &histogram[c0][c1min][c2];
				for (c1 = c1min; c1 <= c1max; c1++, histp += HIST_C2_ELEMS)
					if (*histp != 0) {
						boxp->c2max = c2max = c2;
						goto have_c2max;
					}
			}
have_c2max:

	/* Update box volume.
	* We use 2-norm rather than real volume here; this biases the method
	* against making long narrow boxes, and it has the side benefit that
	* a box is splittable iff norm > 0.
	* Since the differences are expressed in histogram-cell units,
	* we have to shift back to JSAMPLE units to get consistent distances;
	* after which, we scale according to the selected distance scale factors.
	*/
	dist0 = ((c0max - c0min) << C0_SHIFT) * C0_SCALE;
	dist1 = ((c1max - c1min) << C1_SHIFT) * C1_SCALE;
	dist2 = ((c2max - c2min) << C2_SHIFT) * C2_SCALE;
	boxp->volume = dist0*dist0 + dist1*dist1 + dist2*dist2;

	/* Now scan remaining volume of box and compute population */
	ccount = 0;
	for (c0 = c0min; c0 <= c0max; c0++)
		for (c1 = c1min; c1 <= c1max; c1++) {
			histp = &histogram[c0][c1][c2min];
			for (c2 = c2min; c2 <= c2max; c2++, histp++)
				if (*histp != 0) {
					ccount++;
				}
		}
	boxp->colorcount = ccount;
}


int
median_cut(j_decompress_ptr cinfo, boxptr boxlist, int numboxes,
	int desired_colors)
	/* Repeatedly select and split the largest box until we have enough boxes */
{
	int n, lb;
	int c0, c1, c2, cmax;
	boxptr b1, b2;

	while (numboxes < desired_colors) {
		/* Select box to split.
		* Current algorithm: by population for first half, then by volume.
		*/
		if ((numboxes * 2) <= desired_colors) {
			b1 = find_biggest_color_pop(boxlist, numboxes);
		}
		else {
			b1 = find_biggest_volume(boxlist, numboxes);
		}
		if (b1 == NULL)     /* no splittable boxes left! */
			break;
		b2 = &boxlist[numboxes];    /* where new box will go */
									/* Copy the color bounds to the new box. */
		b2->c0max = b1->c0max; b2->c1max = b1->c1max; b2->c2max = b1->c2max;
		b2->c0min = b1->c0min; b2->c1min = b1->c1min; b2->c2min = b1->c2min;
		/* Choose which axis to split the box on.
		* Current algorithm: longest scaled axis.
		* See notes in update_box about scaling distances.
		*/
		c0 = ((b1->c0max - b1->c0min) << C0_SHIFT) * C0_SCALE;
		c1 = ((b1->c1max - b1->c1min) << C1_SHIFT) * C1_SCALE;
		c2 = ((b1->c2max - b1->c2min) << C2_SHIFT) * C2_SCALE;
		/* We want to break any ties in favor of green, then red, blue last.
		* This code does the right thing for R,G,B or B,G,R color orders only.
		*/
#if RGB_RED == 0
		cmax = c1; n = 1;
		if (c0 > cmax) { cmax = c0; n = 0; }
		if (c2 > cmax) { n = 2; }
#else
		cmax = c1; n = 1;
		if (c2 > cmax) { cmax = c2; n = 2; }
		if (c0 > cmax) { n = 0; }
#endif

		/* Choose split point along selected axis, and update box bounds.
		* Current algorithm: split at halfway point.
		* (Since the box has been shrunk to minimum volume,
		* any split will produce two nonempty subboxes.)
		* Note that lb value is max for lower box, so must be < old max.
		*/
		switch (n) {
		case 0:
			lb = (b1->c0max + b1->c0min) / 2;
			b1->c0max = lb;
			b2->c0min = lb + 1;
			break;
		case 1:
			lb = (b1->c1max + b1->c1min) / 2;
			b1->c1max = lb;
			b2->c1min = lb + 1;
			break;
		case 2:
			lb = (b1->c2max + b1->c2min) / 2;
			b1->c2max = lb;
			b2->c2min = lb + 1;
			break;
		}
		/* Update stats for boxes */
		update_box(cinfo, b1);
		update_box(cinfo, b2);
		numboxes++;
	}
	return numboxes;
}


void
compute_color(j_decompress_ptr cinfo, boxptr boxp, int icolor)
/* Compute representative color for a box, put it in colormap[icolor] */
{
	/* Current algorithm: mean weighted by pixels (not colors) */
	/* Note it is important to get the rounding correct! */
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	hist3d histogram = cquantize->histogram;
	histptr histp;
	int c0, c1, c2;
	int c0min, c0max, c1min, c1max, c2min, c2max;
	long count;
	long total = 0;
	long c0total = 0;
	long c1total = 0;
	long c2total = 0;

	c0min = boxp->c0min;  c0max = boxp->c0max;
	c1min = boxp->c1min;  c1max = boxp->c1max;
	c2min = boxp->c2min;  c2max = boxp->c2max;

	for (c0 = c0min; c0 <= c0max; c0++)
		for (c1 = c1min; c1 <= c1max; c1++) {
			histp = &histogram[c0][c1][c2min];
			for (c2 = c2min; c2 <= c2max; c2++) {
				if ((count = *histp++) != 0) {
					total += count;
					c0total += ((c0 << C0_SHIFT) + ((1 << C0_SHIFT) >> 1)) * count;
					c1total += ((c1 << C1_SHIFT) + ((1 << C1_SHIFT) >> 1)) * count;
					c2total += ((c2 << C2_SHIFT) + ((1 << C2_SHIFT) >> 1)) * count;
				}
			}
		}

	cinfo->colormap[0][icolor] = (JSAMPLE)((c0total + (total >> 1)) / total);
	cinfo->colormap[1][icolor] = (JSAMPLE)((c1total + (total >> 1)) / total);
	cinfo->colormap[2][icolor] = (JSAMPLE)((c2total + (total >> 1)) / total);
}

static void
select_colors(j_decompress_ptr cinfo, int desired_colors)
/* Master routine for color selection */
{
	boxptr boxlist;
	int numboxes;
	int i;

	/* Allocate workspace for box list */
	boxlist = (boxptr)malloc(desired_colors * sizeof(box));
	/* Initialize one box containing whole space */
	numboxes = 1;
	boxlist[0].c0min = 0;
	boxlist[0].c0max = MAXJSAMPLE >> C0_SHIFT;
	boxlist[0].c1min = 0;
	boxlist[0].c1max = MAXJSAMPLE >> C1_SHIFT;
	boxlist[0].c2min = 0;
	boxlist[0].c2max = MAXJSAMPLE >> C2_SHIFT;
	/* Shrink it to actually-used volume and set its statistics */
	update_box(cinfo, &boxlist[0]);
	/* Perform median-cut to produce final box list */
	numboxes = median_cut(cinfo, boxlist, numboxes, desired_colors);
	/* Compute the representative color for each box, fill colormap */
	for (i = 0; i < numboxes; i++)
		compute_color(cinfo, &boxlist[i], i);
	cinfo->actual_number_of_colors = numboxes;

	free(boxlist); //FIXME?? I don't know if this is correct - VS
}



/*
* Finish up at the end of each pass.
*/

void
finish_pass1(j_decompress_ptr cinfo)
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;

	/* Select the representative colors and fill in cinfo->colormap */
	cinfo->colormap = cquantize->sv_colormap;
	select_colors(cinfo, cquantize->desired);
	/* Force next pass to zero the color index table */
	cquantize->needs_zeroed = true;
}

void
finish_pass2(j_decompress_ptr WXUNUSED(cinfo))
{
	/* no work */
}
void
prescan_quantize(j_decompress_ptr cinfo, JSAMPARRAY input_buf,
	JSAMPARRAY WXUNUSED(output_buf), int num_rows)
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	JSAMPROW ptr;
	histptr histp;
	hist3d histogram = cquantize->histogram;
	int row;
	JDIMENSION col;
	JDIMENSION width = cinfo->output_width;

	for (row = 0; row < num_rows; row++) {
		ptr = input_buf[row];
		for (col = width; col > 0; col--) {

			{

				/* get pixel value and index into the histogram */
				histp = &histogram[GETJSAMPLE(ptr[0]) >> C0_SHIFT]
					[GETJSAMPLE(ptr[1]) >> C1_SHIFT]
				[GETJSAMPLE(ptr[2]) >> C2_SHIFT];
				/* increment, check for overflow and undo increment if so. */
				if (++(*histp) <= 0)
					(*histp)--;
			}
			ptr += 3;
		}
	}
}

void
start_pass_2_quant(j_decompress_ptr cinfo, bool is_pre_scan)
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;
	hist3d histogram = cquantize->histogram;

	if (is_pre_scan) {
		/* Set up method pointers */
		cquantize->pub.color_quantize = prescan_quantize;
		cquantize->pub.finish_pass = finish_pass1;
		cquantize->needs_zeroed = true; /* Always zero histogram */
	}
	else {
		/* Set up method pointers */
		cquantize->pub.color_quantize = pass2_fs_dither;
		cquantize->pub.finish_pass = finish_pass2;

		{
			size_t arraysize = (size_t)((cinfo->output_width + 2) *
				(3 * sizeof(FSERROR)));
			/* Allocate Floyd-Steinberg workspace if we didn't already. */
			if (cquantize->fserrors == NULL)
				cquantize->fserrors = (wxInt16*)malloc(arraysize);
			/* Initialize the propagated errors to zero. */
			memset((void  *)cquantize->fserrors, 0, arraysize);
			/* Make the error-limit table if we didn't already. */
			if (cquantize->error_limiter == NULL)
				init_error_limit(cinfo);
			cquantize->on_odd_row = false;
		}

	}
	/* Zero the histogram or inverse color map, if necessary */
	if (cquantize->needs_zeroed) {
		for (int i = 0; i < HIST_C0_ELEMS; i++) {
			memset((void  *)histogram[i], 0,
				HIST_C1_ELEMS*HIST_C2_ELEMS * sizeof(histcell));
		}
		cquantize->needs_zeroed = false;
	}
}


void
new_color_map_2_quant(j_decompress_ptr cinfo)
{
	my_cquantize_ptr cquantize = (my_cquantize_ptr)cinfo->cquantize;

	/* Reset the inverse color map */
	cquantize->needs_zeroed = true;
}

void
jinit_2pass_quantizer(j_decompress_ptr cinfo)
{
	my_cquantize_ptr cquantize;
	int i;

	cquantize = (my_cquantize_ptr)malloc(sizeof(my_cquantizer));
	cinfo->cquantize = (jpeg_color_quantizer *)cquantize;
	cquantize->pub.start_pass = start_pass_2_quant;
	cquantize->pub.new_color_map = new_color_map_2_quant;
	cquantize->fserrors = NULL;   /* flag optional arrays not allocated */
	cquantize->error_limiter = NULL;


	/* Allocate the histogram/inverse colormap storage */
	cquantize->histogram = (hist3d)malloc(HIST_C0_ELEMS * sizeof(hist2d));
	for (i = 0; i < HIST_C0_ELEMS; i++) {
		cquantize->histogram[i] = (hist2d)malloc(HIST_C1_ELEMS*HIST_C2_ELEMS * sizeof(histcell));
	}
	cquantize->needs_zeroed = true; /* histogram is garbage now */

									/* Allocate storage for the completed colormap, if required.
									* We do this now since it is  storage and may affect
									* the memory manager's space calculations.
									*/
	{
		/* Make sure color count is acceptable */
		int desired = cinfo->desired_number_of_colors;

		cquantize->sv_colormap = (JSAMPARRAY)malloc(sizeof(JSAMPROW) * 3);
		cquantize->sv_colormap[0] = (JSAMPROW)malloc(sizeof(JSAMPLE) * desired);
		cquantize->sv_colormap[1] = (JSAMPROW)malloc(sizeof(JSAMPLE) * desired);
		cquantize->sv_colormap[2] = (JSAMPROW)malloc(sizeof(JSAMPLE) * desired);

		cquantize->desired = desired;
	}

	/* Allocate Floyd-Steinberg workspace if necessary.
	* This isn't really needed until pass 2, but again it is  storage.
	* Although we will cope with a later change in dither_mode,
	* we do not promise to honor max_memory_to_use if dither_mode changes.
	*/
	{
		cquantize->fserrors = (FSERRPTR)malloc(
			(size_t)((cinfo->output_width + 2) * (3 * sizeof(FSERROR))));
		/* Might as well create the error-limiting table too. */
		init_error_limit(cinfo);
	}
}

void
prepare_range_limit_table(j_decompress_ptr cinfo)
/* Allocate and fill in the sample_range_limit table */
{
	JSAMPLE * table;
	int i;

	table = (JSAMPLE *)malloc((5 * (MAXJSAMPLE + 1) + CENTERJSAMPLE) * sizeof(JSAMPLE));
	cinfo->srl_orig = table;
	table += (MAXJSAMPLE + 1);  /* allow negative subscripts of simple table */
	cinfo->sample_range_limit = table;
	/* First segment of "simple" table: limit[x] = 0 for x < 0 */
	memset(table - (MAXJSAMPLE + 1), 0, (MAXJSAMPLE + 1) * sizeof(JSAMPLE));
	/* Main part of "simple" table: limit[x] = x */
	for (i = 0; i <= MAXJSAMPLE; i++)
		table[i] = (JSAMPLE)i;
	table += CENTERJSAMPLE;   /* Point to where post-IDCT table starts */
							  /* End of simple table, rest of first half of post-IDCT table */
	for (i = CENTERJSAMPLE; i < 2 * (MAXJSAMPLE + 1); i++)
		table[i] = MAXJSAMPLE;
	/* Second half of post-IDCT table */
	memset(table + (2 * (MAXJSAMPLE + 1)), 0,
		(2 * (MAXJSAMPLE + 1) - CENTERJSAMPLE) * sizeof(JSAMPLE));
	memcpy(table + (4 * (MAXJSAMPLE + 1) - CENTERJSAMPLE),
		cinfo->sample_range_limit, CENTERJSAMPLE * sizeof(JSAMPLE));
}

void DoQuantize(unsigned w, unsigned h, unsigned char **in_rows, unsigned char **out_rows,
	unsigned char *palette, int desiredNoColours)
{
	j_decompress dec;
	my_cquantize_ptr cquantize;

	dec.output_width = w;
	dec.desired_number_of_colors = desiredNoColours;
	prepare_range_limit_table(&dec);
	jinit_2pass_quantizer(&dec);
	cquantize = (my_cquantize_ptr)dec.cquantize;


	cquantize->pub.start_pass(&dec, true);
	cquantize->pub.color_quantize(&dec, in_rows, out_rows, h);
	cquantize->pub.finish_pass(&dec);

	cquantize->pub.start_pass(&dec, false);
	cquantize->pub.color_quantize(&dec, in_rows, out_rows, h);
	cquantize->pub.finish_pass(&dec);


	for (int i = 0; i < dec.desired_number_of_colors; i++) {
		palette[3 * i + 0] = dec.colormap[0][i];
		palette[3 * i + 1] = dec.colormap[1][i];
		palette[3 * i + 2] = dec.colormap[2][i];
	}

	for (int ii = 0; ii < HIST_C0_ELEMS; ii++) free(cquantize->histogram[ii]);
	free(cquantize->histogram);
	free(dec.colormap[0]);
	free(dec.colormap[1]);
	free(dec.colormap[2]);
	free(dec.colormap);
	free(dec.srl_orig);

	//free(cquantize->error_limiter);
	free((void*)(cquantize->error_limiter - MAXJSAMPLE)); // To reverse what was done to it

	free(cquantize->fserrors);
	free(cquantize);
}
