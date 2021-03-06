/* *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* */
/* ** Copyright UCAR (c) 1992 - 2001 */
/* ** University Corporation for Atmospheric Research(UCAR) */
/* ** National Center for Atmospheric Research(NCAR) */
/* ** Research Applications Program(RAP) */
/* ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA */
/* ** 2001/11/19 23:5:23 */
/* *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* */
#ifndef SWAP_H
#define SWAP_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Module: swap.h
 *
 * Author: Mike Dixon
 *	
 * Date: Jan 1996
 *
 */

/*
 * Swap module
 *
 * Description:
 *
 * The routines in this module swap arrays of 16-bit, 32-bit and 64-bit
 * words.
 *
 */

#include <dataport/port_types.h>

/**********************************************************************
 * SWAP_array_64()
 *
 * Performs an in-place 64-bit value byte swap.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.
 *
 */

extern ui32 SWAP_array_64(void *array, ui32 nbytes);

/**********************************************************************
 * SWAP_array_32()
 *
 * Performs an in-place 32-bit word byte swap.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.
 *
 */

extern ui32 SWAP_array_32(void *array, ui32 nbytes);

/**********************************************************************
 * SWAP_array_16()
 *
 * Performs an in-place 16-bit word byte swap.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.
 *
 */

extern ui32 SWAP_array_16(void *array, ui32 nbytes);

/* See dataport/port_types.h for whish OS'es to include */
/* si64 and ui64 not defined for the others (yet) */
#if defined(IRIX6) || defined(DECOSF1) || defined(SUNOS5_64) 
/********************************
 *  SWAP_si64
 *  Swaps a single si64
 */

extern si64 SWAP_si64(si64 x);

/********************************
 *  SWAP_ui64
 *  Swaps a single ui64
 */

extern ui64 SWAP_ui64(ui64 x);
#endif

/********************************
 *  SWAP_fl64
 *  Swaps a single fl64
 */

extern fl64 SWAP_fl64(fl64 x);

/********************************
 *  SWAP_si32
 *  Swaps a single si32
 */

extern si32 SWAP_si32(si32 x);

/********************************
 *  SWAP_ui32
 *  Swaps a single ui32
 */

extern ui32 SWAP_ui32(ui32 x);

/********************************
 *  SWAP_fl32
 *  Swaps a single fl32
 */

extern fl32 SWAP_fl32(fl32 x);

/********************************
 *  SWAP_si16
 *  Swaps a single si16
 */

extern si16 SWAP_si16(si16 x);

/********************************
 *  SWAP_ui16
 *  Swaps a single ui16
 */

extern ui16 SWAP_ui16(ui16 x);

#ifdef __cplusplus
}
#endif

#endif /* SWAP_H */

