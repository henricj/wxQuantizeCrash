#include "stdafx.h"
#include "quantize.h"

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
