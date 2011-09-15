/*
 *  crlibm_private.h
 *  
 * This file contains useful tools and data for the crlibm functions.
 *
 */

#ifndef CRLIBM_PRIVATE_H
#define CRLIBM_PRIVATE_H 1

#define SCS_NB_WORDS 8

#include <inttypes.h>

typedef union {
   int32_t i[2];
   int64_t l;
   double d;
} db_number;

struct scs {
  /** the digits, as 32 bits words */
    uint32_t h_word[SCS_NB_WORDS];
    /** Used to store Nan,+/-0, Inf, etc and then let the hardware handle them */
    db_number exception;
          /** This corresponds to the exponent in an FP format, but here we are 
        in base 2^32  */
   int index;
    /** The sign equals 1 or -1*/
     int sign;
};


typedef struct scs scs;

#define HI 1
#define LO 0

#define AVOID_BRANCHES 1


/*
 * In the following, when an operator is preceded by a '@' it means that we
 * are considering the IEEE-compliant machine operator, otherwise it
 * is the mathematical operator.
 *
 */

#define Add12Cond(s, r, a, b)      \
{                                  \
    double _u1, _u2, _u3, _u4;     \
    double  _a=a, _b=b;            \
                                   \
    s = _a + _b;                   \
    _u1 = s - _a;                  \
    _u2 = s - _u1;                 \
    _u3 = _b - _u1;                \
    _u4 = _a - _u2;                \
    r = _u4 + _u3;                 \
} 


/*
 *  computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 * under the condition  a >= b
 */
#define Add12(s, r, a, b)          \
        {double _z, _a=a, _b=b;    \
         s = _a + _b;              \
         _z = s - _a;              \
         r = _b - _z;   }            

/*
 * Functions to computes double-double addition: zh+zl = xh+xl + yh+yl
 * knowing that xh>yh
 * relative error is smaller than 2^-103 
 */


#if AVOID_BRANCHES
#define Add22Cond(zh,zl,xh,xl,yh,yl)                                                   \
do {                                                                                   \
  double _v1, _v2, _v3, _v4;                                                           \
                                                                                       \
  Add12Cond(_v1, _v2, (xh), (yh));                                                     \
  _v3 = (xl) + (yl);                                                                   \
  _v4 = _v2 + _v3;                                                                     \
  Add12((*(zh)),(*(zl)),_v1,_v4);                                                      \
} while (2+2==5) 
#else
#define Add22Cond(zh,zl,xh,xl,yh,yl)                                                   \
do {                                                                                   \
  double _r,_s;                                                                        \
  _r = (xh)+(yh);                                                                      \
  _s = ((ABS(xh)) > (ABS(yh)))? ((xh)-_r+(yh)+(yl)+(xl)) : ((yh)-_r+(xh)+(xl)+(yl));   \
  *zh = _r+_s;                                                                         \
  *zl = (_r - (*zh)) + _s;                                                             \
} while(2+2==5)
#endif
  

#define Add22(zh,zl,xh,xl,yh,yl)         \
do {                                     \
double _r,_s;                            \
_r = (xh)+(yh);                          \
_s = ((((xh)-_r) +(yh)) + (yl)) + (xl);  \
*zh = _r+_s;                             \
*zl = (_r - (*zh)) + _s;                 \
} while(0)


/*
 * computes rh and rl such that rh + rl = a * b with rh = a @* b exactly
 * under the conditions : a < 2^970 et b < 2^970 
 */
#define Mul12(rh,rl,u,v)                        \
{                                               \
  const double c  = 134217729.; /* 2^27 +1 */   \
  double up, u1, u2, vp, v1, v2;                \
  double _u=u, _v=v;                            \
  up = _u*c;        vp = _v*c;                  \
  u1 = (_u-up)+up;  v1 = (_v-vp)+vp;            \
  u2 = _u-u1;       v2 = _v-v1;                 \
                                                \
  *rh = _u*_v;                                  \
  *rl = (((u1*v1-*rh)+(u1*v2))+(u2*v1))+(u2*v2);\
}

/*
  double _u =u, _v=v;                           \
 __m128d _u_v = _mm_set_pd(_u, _v);            \
*/                                                \
/*
 * Computes rh and rl such that rh + rl = a * b and rh = a @* b exactly
 */
#define Mul12Cond(rh, rl, a,  b) \
{\
  const double two_em53 = 1.1102230246251565404e-16; /* 0x3CA00000, 0x00000000 */\
  const double two_e53  = 9007199254740992.;         /* 0x43400000, 0x00000000 */\
  double u, v;                                               \
  db_number _a=a, _b=b;                                      \
                                                             \
  if (_a.i[HI]>0x7C900000) u = _a*two_em53;                  \
  else            u = _a;                                    \
  if (_b.i[HI]>0x7C900000) v = _b*two_em53;                  \
  else            v = _b;                                    \
                                                             \
  Mul12(rh, rl, u, v);                                       \
                                                             \
  if (_a.i[HI]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
  if (_b.i[HI]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
}



/*
 * computes double-double multiplication: zh+zl = (xh+xl) *  (yh+yl)
 * relative error is smaller than 2^-102
 */
  

  
#define Mul22(zh,zl,xh,xl,yh,yl)                      \
{                                                     \
double mh, ml;                                        \
						      \
  const double c = 134217729.;			      \
  double up, u1, u2, vp, v1, v2;		      \
						      \
  up = (xh)*c;        vp = (yh)*c;		      \
  u1 = ((xh)-up)+up;  v1 = ((yh)-vp)+vp;	      \
  u2 = (xh)-u1;       v2 = (yh)-v1;                   \
  						      \
  mh = (xh)*(yh);				      \
  ml = (((u1*v1-mh)+(u1*v2))+(u2*v1))+(u2*v2);	      \
						      \
  ml += (xh)*(yl) + (xl)*(yh);			      \
  *zh = mh+ml;					      \
  *zl = mh - (*zh) + ml;                              \
}


/* Additional double-double operators */

/* Eps Mul122 <= 2^-102 */
#define Mul122(resh,resl,a,bh,bl)                 \
{                                                 \
    double _t1, _t2, _t3, _t4;                    \
                                                  \
    Mul12(&_t1,&_t2,(a),(bh));                    \
    _t3 = (a) * (bl);                             \
    _t4 = _t2 + _t3;                              \
    Add12((*(resh)),(*(resl)),_t1,_t4);           \
}

/* Eps MulAdd212 <= 2^-100 for |a * (bh + bl)| <= 1/4 * |ch + cl| */
#define MulAdd212(resh,resl,ch,cl,a,bh,bl)           \
{                                                    \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8;   \
                                                     \
    Mul12(&_t1,&_t2,(a),(bh));                       \
    Add12(_t3,_t4,(ch),_t1);                         \
    _t5 = (bl) * (a);                                \
    _t6 = (cl) + _t2;                                \
    _t7 = _t5 + _t6;                                 \
    _t8 = _t7 + _t4;                                 \
    Add12((*(resh)),(*(resl)),_t3,_t8);              \
}

/* Eps MulAdd212 <= 2^-100 
   for |(ah + bh) * (bh + bl)| <= 1/4 * |ch + cl| 
*/
#define MulAdd22(resh,resl,ch,cl,ah,al,bh,bl)        \
{                                                    \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8;   \
    double _t9, _t10;                                \
                                                     \
    Mul12(&_t1,&_t2,(ah),(bh));                      \
    Add12(_t3,_t4,(ch),_t1);                         \
    _t5 = (ah) * (bl);                               \
    _t6 = (al) * (bh);                               \
    _t7 = _t2 + (cl);                                \
    _t8 = _t4 + _t7;                                 \
    _t9 = _t5 + _t6;                                 \
    _t10 = _t8 + _t9;                                \
    Add12((*(resh)),(*(resl)),_t3,_t10);             \
}

#define Add122(resh,resl,a,bh,bl)                    \
{                                                    \
    double _t1, _t2, _t3;                            \
                                                     \
    Add12(_t1,_t2,(a),(bh));                         \
    _t3 = _t2 + (bl);                                \
    Add12((*(resh)),(*(resl)),_t1,_t3);              \
}    

#define Add122Cond(resh,resl,a,bh,bl)                \
{                                                    \
    double _t1, _t2, _t3;                            \
                                                     \
    Add12Cond(_t1,_t2,(a),(bh));                     \
    _t3 = _t2 + (bl);                                \
    Add12((*(resh)),(*(resl)),_t1,_t3);              \
}    


#define Add212(resh,resl,ah,al,b)                    \
{                                                    \
    double _t1, _t2, _t3;                            \
                                                     \
    Add12(_t1,_t2,(ah),b);                           \
    _t3 = _t2 + (al);                                \
    Add12((*(resh)),(*(resl)),_t1,_t3);              \
}


/* In the following the one-line computation of _cl was split so that
   icc(8.1) would compile it properly. It's a bug of icc */

#if DEKKER_AS_FUNCTIONS
extern void Div22(double *z, double *zz, double x, double xx, double y, double yy);
#else
#define  Div22(pzh,pzl,xh,xl,yh,yl)  {           \
  double _ch,_cl,_uh,_ul;                        \
  _ch=(xh)/(yh);   Mul12(&_uh,&_ul,_ch,(yh));    \
  _cl=((xh)-_uh);                                \
  _cl -= _ul;                                    \
  _cl += (xl);                                   \
  _cl -= _ch*(yl);                               \
  _cl /= (yh);                                   \
  *pzh=_ch+_cl;   *pzl=(_ch-(*pzh))+_cl;         \
}
#endif /* DEKKER_AS_FUNCTIONS */

#endif /*CRLIBM_PRIVATE_H*/
