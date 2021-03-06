/* *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* */
/* ** Copyright UCAR (c) 1992 - 2001 */
/* ** University Corporation for Atmospheric Research(UCAR) */
/* ** National Center for Atmospheric Research(NCAR) */
/* ** Research Applications Program(RAP) */
/* ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA */
/* ** 2001/11/19 23:5:23 */
/* *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* */
#ifndef SMALLEND_H
#define SMALLEND_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Module: smallend.h
 *
 * Author: Gerry Wiener
 *	
 * Date:   Jan 26, 1995
 *
 * Modified by Mike Dixon, Jan 1996.
 *
 */

/*
 * Smallend module
 *
 * Description:
 *   This module consists of library routines to convert integers and
 * floating point values to and from small endian format.  Small endian
 * format assigns lower order bytes to storage with smaller addresses.
 * For example, if the number 1 is to be stored using a two byte small
 * endian integer, it would be represented as
 *
 * address n+1    address n
 * 00000000       00000001
 *
 * The same number would be stored as
 *
 * address n+1    address n
 * 00000001       00000000
 *
 * in big endian format.  In the following code, we assume that the
 * underlying machine uses either big endian or little endian addressing.
 * If this is not the case, the module has to be rewritten to support
 * the underlying addressing scheme.
 *   The strategy of this module is to provide tools for machine
 * independent byte storage.  The routines in this module were designed
 * for efficiency, utility and portability .  The routines do not pad but
 * assume that different types will be converted to appropriate fixed
 * storage sizes.  In order to use these routines, one needs to determine
 * which C integer types on the underlying machine have 8 bits, 16 bits
 * and 32 bits.  On current machines one would use unsigned char,
 * unsigned short and unsigned int.  One would then reset ui08,
 * ui16, ui32, in smallend.h appropriately.
 *
 * Note that as of this date, Jan 24, 1995, Unix workstations such as
 * Sun, SGI, Digital, IBM, HP implement the following sizes:
 *
 * char   -> 1 byte
 * short  -> 2 bytes
 * int    -> 4 bytes
 * long   -> 4 bytes or 8 bytes (Digital)
 * float  -> 4 bytes
 * double -> 8 bytes
 *
 * MSDOS/WINDOWS typically assign 2 bytes to integers and 4 bytes to longs.
 *
 * IMPORTANT NOTE: The software assumes that the floating point
 * implementation is identical on all machines in question except for
 * byte ordering.
 *
 * In the future, it may be the case that types are assigned larger byte
 * lengths perhaps 12 bytes or 16 bytes.  In such cases, these routines
 * could be extended.
 */

#include <dataport/port_types.h>

/*
 * determine whether the underlying machine is small endian or not
 */

extern int SE_is_small_endian(void);

/*
 * SE_reverse()
 *
 * Reverses the sense of this library. Therefore,
 * is called once, BigEndian values are set.
 * If called twice, goes back to SmallEndian.
 */

extern void SE_reverse(void);

/*
 * in-place array converting
 */

/****************************************************
 * SE_swap_array_64()
 *
 * Performs an in-place 64-bit value byte swap.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.
 */

extern si32 SE_swap_array_64(void *array, ui32 nbytes);

/**********************************************************************
 * SE_swap_array_32()
 *
 * Performs an in-place 32-bit word byte swap, if necessary, to produce
 * SE representation from machine representation, or vice-versa.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 *
 */

extern si32 SE_swap_array_32(void *array, ui32 nbytes);

/**********************************************************************
 * SE_swap_array_16()
 *
 * Performs an in-place 16-bit word byte swap, if necessary, to produce
 * SE representation from machine representation, or vice-versa.
 *
 * Array must be aligned.
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 *
 */

extern si32 SE_swap_array_16(void *array, ui32 nbytes);

/********************
 * SE_from_array_32()
 * Converts an array of 32's
 * Replaces XDRU_fromhl
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_from_array_32(void *array, ui32 nbytes);
     
/******************
 * SE_to_array_32()
 * Converts an array of 32's
 * Replaces XDRU_tohl
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_to_array_32(void *array, ui32 nbytes);
     
/********************
 * SE_from_array_64()
 * Converts an array of 64's
 * Replaces XDRU_fromhl
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_from_array_64(void *array, ui32 nbytes);
     
/******************
 * SE_to_array_64()
 * Converts an array of 64's
 * Replaces XDRU_tohl
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_to_array_64(void *array, ui32 nbytes);
     
/********************
 * SE_from_array_16()
 * Converts an array of 16's
 * Replaces XDRU_fromhs
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_from_array_16(void *array, ui32 nbytes);

/******************
 * SE_to_array_16()
 * Converts an array of 16's
 * Replaces XDRU_tohs
 *
 * Returns the number of bytes converted.  Note that this number is 0
 * if the local machine is small-endian.
 */

extern si32 SE_to_array_16(void *array, ui32 nbytes);

/********************************
 *  SE_from_si32 replaces htonl()
 *  Converts a single si32
 *
 *  Returns the converted number.
 */

extern si32 SE_from_si32(si32 x);

/******************************
 *  SE_to_si32 replaces ntohl()
 *  Converts a single si32
 *
 *  Returns the converted number.
 */

extern si32 SE_to_si32(si32 x);

/********************************
 *  SE_from_si16 replaces htons()
 *  Converts a single si16
 *
 *  Returns the converted number.
 */

extern si16 SE_from_si16(si16 x);

/******************************
 *  SE_to_si16 replaces ntohs()
 *  Converts a single si16
 *
 *  Returns the converted number.
 */

extern si16 SE_to_si16(si16 x);

/********************************
 *  SE_from_ui32
 *  Converts a single ui32
 *
 *  Returns the converted number.
 */

extern ui32 SE_from_ui32(ui32 x);

/******************************
 *  SE_to_ui32
 *  Converts a single ui32
 *
 *  Returns the converted number.
 */

extern ui32 SE_to_ui32(ui32 x);

/********************************
 *  SE_from_fl32
 *  Converts a single fl32
 *
 *  Returns the converted number.
 */

extern fl32 SE_from_fl32(fl32 x);

/******************************
 *  SE_to_fl32
 *  Converts a single fl32
 *
 *  Returns the converted number.
 */

extern fl32 SE_to_fl32(fl32 x);

/********************************
 *  SE_from_ui16
 *  Converts a single ui16
 *
 *  Returns the converted number.
 */

extern ui16 SE_from_ui16(ui16 x);

/******************************
 *  SE_to_ui16
 *  Converts a single ui16
 *
 *  Returns the converted number.
 */

extern ui16 SE_to_ui16(ui16 x);

#ifdef __cplusplus
}
#endif

#endif /* SMALLEND_H */

