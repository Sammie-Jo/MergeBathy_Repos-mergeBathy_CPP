#ifndef __UTILITY_NVDEFS_H__
#define __UTILITY_NVDEFS_H__

#ifdef  __cplusplus
extern "C" {
#endif


/* WGS-84 semi-major & semi-minor axes, resp. */

#define         NV_A0           6378137.0L
#define         NV_B0           6356752.314245L

#define         NV_RF           298.257223563L

#define         LOG2            0.3010299956639812L

#define         NV_DEG_TO_RAD   0.017453293L
#define         NV_RAD_TO_DEG   57.2957795147195L


#define         NINT(a)         ((a) < 0.0 ? (int) ((a) - 0.5) : (int) ((a) + 0.5))


#ifndef SIGN
  #define       SIGN(a)         ((a) < 0.0 ? -1.0 : 1.0)
#endif


#ifndef MAX
  #define       MAX(x,y)        (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
  #define       MIN(x,y)        (((x) < (y)) ? (x) : (y))
#endif


#define         DPRINT          fprintf (stderr, "%s %d\n", __FILE__, __LINE__);fflush (stderr);



#ifdef  __cplusplus
}
#endif

#endif
