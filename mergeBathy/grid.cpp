/**********************************************************************
* CC0 License
**********************************************************************
* MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
* Modified in 2015 by Samantha J.Zambo(samantha.zambo@gmail.com) while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Todd Holland while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Nathaniel Plant while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Kevin Duvieilh while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Paul Elmore while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Will Avera while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Brian Bourgeois while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by A.Louise Perkins while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by David Lalejini while employed by the U.S.Naval Research Laboratory.
* To the extent possible under law, the author(s) and the U.S.Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide.This software is distributed without any warranty.
* You should have received a copy of the CC0 Public Domain Dedication along with this software.If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
**********************************************************************/
//============================================================================
// Grid.C
// Downloaded from http:\\oldmill.uchicago.edu/~wilder/Code/grid/ << broken
// Working site as of 3 Mar 2018: https://sites.google.com/site/jivsoft/Home

#include "grid.h"

using namespace std;

//============================================================================
// SPECIALIZATIONS FOR UCGRID

template<> ostream& operator<<(ostream& os, const ucgrid& g)
{
  for (uint i = 0; i < g.rows(); ++i) {
    for (uint j = 0; j < g.cols(); ++j)
      os << setw(4) << (int) g(i, j) << " ";
    os << "\n";
  }
  return os;
}

//============================================================================
// SPECIALIZATIONS FOR DGRID

#ifdef USE_LAPACK

const int izero = 0, ione = 1;
const double dzero = 0, done = 1;
const cplx zzero = 0.0, zone = 1.0;
const char notrans = 'N';

extern "C" {
//============================================================================
// Blas 1 routine to compute y = a * x + y

// Note: dcopy/zcopy and dscal/zscal are not included here as they are
// typically not much, if any, faster than non-Blas routines.

void daxpy_(const int*, const double*, const double*, const int*,
            double*, const int*);
void zaxpy_(const int*, const cplx*, const cplx*, const int*,
            cplx*, const int*);

void drot_(const int*, double*, const int*, double*,
           const int*, const double*, const double*);
void zrot_(const int*, cplx*, const int*, cplx*,
           const int*, const double*, const cplx*);

void dlartg_(const double*, const double*, double*, double*, double*);
void zlartg_(const cplx*, const cplx*, double*, cplx*, cplx*);

//============================================================================
// Blas 3 routines to compute the matrix multiplication product of matrices.

void dtrmm_(const char* side, const char* uplo, const char* tra,
            const char* diag, const int* m, const int* n, const double* alpha,
	    const double* a, const int* lda, double* b, const int* ldb);
void ztrmm_(const char* side, const char* uplo, const char* tra,
            const char* diag, const int* m, const int* n, const cplx* alpha,
	    const cplx* a, const int* lda, cplx* b, const int* ldb);
void dgemm_(const char*, const char*, const int*, const int*, const int*,
            const double*, const double*, const int*, const double*,
            const int*, const double*, double*, const int*);
void zgemm_(const char*, const char*, const int*, const int*, const int*,
            const cplx*, const cplx*, const int*, const cplx*,
            const int*, const cplx*, cplx*, const int*);

//============================================================================
// Blas 3 routines to solve a AX=B where A is triangular.
// These are similar to the Lapack trtrs routines, but they also
// allow for solving XA=B.

void dtrsm_(const char*, const char*, const char*, const char*,
            const int*, const int*, const double*, const double*,
            const int*, const double*, const int*);
void ztrsm_(const char*, const char*, const char*, const char*,
            const int*, const int*, const cplx*, const cplx*,
            const int*, const cplx*, const int*);

//============================================================================
// Lapack routines to solve overdetermined or underdetermined systems.
// using a QR or LQ factorization.

void dgels_(const char* trans, const int* m, const int* n, const int* nrhs,
            double* a, const int* lda, double* b, const int* ldb,
            double* work, const int* lwork, int* info);
void zgels_(const char* trans, const int* m, const int* n, const int* nrhs,
            cplx* a, const int* lda, cplx* b, const int* ldb,
            cplx* work, const int* lwork, int* info);

//============================================================================
// Lapack routines to solve AX=B for square matrices A.

void dgesv_(const int* n, const int* nrhs, double* a, const int* lda,
            int* ipiv, double* b, const int* ldb, int* info);
void zgesv_(const int* n, const int* nrhs, cplx* a, const int* lda,
            int* ipiv, cplx* b, const int* ldb, int* info);
void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a,
            const int* lda, int* ipiv, double* b, const int* ldb,
            double* work, const int* lwork, int* info);
void zsysv_(const char* uplo, const int* n, const int* nrhs, cplx* a,
            const int* lda, int* ipiv, cplx* b, const int* ldb,
            cplx* work, const int* lwork, int* info);
void zhesv_(const char* uplo, const int* n, const int* nrhs, cplx* a,
            const int* lda, int* ipiv, cplx* b, const int* ldb,
            cplx* work, const int* lwork, int* info);
void dposv_(const char* uplo, const int* n, const int* nrhs, double* a,
            const int* lda, double* b, const int* ldb, int* info);
void zposv_(const char* uplo, const int* n, const int* nrhs, cplx* a,
            const int* lda, cplx* b, const int* ldb, int* info);

//============================================================================
// Lapack routines to compute the eigenvalues and (optionally) the
// eigenvectors of a matrix.

void dgeev_(const char* jobvl, const char* jobvr, const int* n, double* a,
            const int* lda, double* wr, double* wi, double* vl,
            const int* ldvl, double* vr, const int* ldvr, double* work,
            const int* lwork, int* info);
void zgeev_(const char* jobvl, const char* jobvr, const int* n, cplx* a,
            const int* lda, cplx* w, cplx* vl, const int* ldvl,
            cplx* vr, const int* ldvr, cplx* work,
            const int* lwork, double* rwork, int* info);
void dsyev_(const char* jobz, const char* uplo, const int* n, double* a,
            const int* lda, double* w, double* work, const int* lwork,
            int* info);
void zheev_(const char* jobz, const char* uplo, const int* n, void* a,
            const int* lda, void* w, void* work, const int* lwork,
            const void* work2, int* info);

//============================================================================
// Lapack routines to compute the singular values and singular vectors
// of a general matrix.

void dgesvd(const char* jobu, const char* jobvt, const int* m, const int* n,
            double* a, const int* lda, double* s, double* u, const int* ldu,
            double* vt, const int* ldvt, double* work, const int* lwork,
            int* info);
void zgesvd(const char* jobu, const char* jobvt, const int* m, const int* n,
            cplx* a, const int* lda, double* s, cplx* u, const int* ldu,
            cplx* vt, const int* ldvt, cplx* work, const int* lwork,
            int* info);

//============================================================================
// Lapack routines to compute the LU/Cholesky Factorization of a matrix
// None are needed for triangular matrices...

void dgetrf_(const int* m, const int* n, double* a, const int* lda,
             int* ipiv, int* info);
void dsytrf_(const char* uplo, const int* n, double* a, const int* lda,
             int* ipiv, double* work, const int* lwork, int* info);
void dpotrf_(const char* uplo, const int* n, double* a,
             const int* lda, int* info);
void zgetrf_(const int* m, const int* n, cplx* a, const int* lda,
             int* ipiv, int* info);
void zsytrf_(const char* uplo, const int* n, cplx* a, const int* lda,
             int* ipiv, cplx* work, const int* lwork, int* info);
void zpotrf_(const char* uplo, const int* n, cplx* a,
             const int* lda, int* info);

//============================================================================
// Lapack routines to solve a linear system using the
// LU/Cholesky Factorization of a matrix

void dgetrs_(const char* trans, const int* n, const int* nrhs,
             const double* a, const int* lda, const int* ipiv,
             double* b, const int* ldb, int* info);
void dsytrs_(const char* uplo, const int* n, const int* nrhs,
             const double* a, const int* lda, const int* ipiv,
             double* b, const int* ldb, int* info);
void dpotrs_(const char* uplo, const int* n, const int* nrhs,
             const double* a, const int* lda, double* b, const int* ldb,
             int* info);
void dtrtrs_(const char* uplo, const char* trans, const char* diag,
             const int* n, const int* nrhs, const double* a,
            const int* lda, double* b, const int* ldb, int* info);
void zgetrs_(const char* trans, const int* n, const int* nrhs,
             const cplx* a, const int* lda, const int* ipiv,
             cplx* b, const int* ldb, int* info);
void zsytrs_(const char* uplo, const int* n, const int* nrhs,
             const cplx* a, const int* lda, const int* ipiv,
             cplx* b, const int* ldb, int* info);
void zpotrs_(const char* uplo, const int* n, const int* nrhs,
             const cplx* a, const int* lda, cplx* b, const int* ldb,
             int* info);
void ztrtrs_(const char* uplo, const char* trans, const char* diag,
             const int* n, const int* nrhs, const cplx* a,
            const int* lda, cplx* b, const int* ldb, int* info);

//============================================================================
// Lapack routines to compute the inverse of a matrix

void dgetri_(const int* n, double* a, const int* lda, const int* ipiv,
             double* work, const int* lwork, int* info);
void dsytri_(const char* uplo, const int* n, double* a,
             const int* lda, const int* ipiv, double* work, int* info);
void dpotri_(const char* uplo, const int* n, double* a,
             const int* lda, int* info);
void dtrtri_(const char* uplo, const char* diag, const int* n,
             double* a, const int* lda, int* info);
void zgetri_(const int* n, cplx* a, const int* lda, const int* ipiv,
             cplx* work, const int* lwork, int* info);
void zsytri_(const char* uplo, const int* n, cplx* a,
             const int* lda, const int* ipiv, cplx* work, int* info);
void zpotri_(const char* uplo, const int* n, cplx* a,
             const int* lda, int* info);
void ztrtri_(const char* uplo, const char* diag, const int* n,
             cplx* a, const int* lda, int* info);
}

//============================================================================
// Compute the matrix product c = a %*% b of two grids/vectors:

// More efficient: matmult(a, b, c); Less efficient: c = matmult(a, b);
// If &c == &a or &c = &b, then 'a' or 'b' will get overwritten.

template<> dgrid&
matmult(const dgrid& a, const dgrid& b, dgrid& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw dgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || &c == &b)
    { dgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    dgemm_(&ta, &tb, &m, &n, &k, &done, &a[0], &ar,
           &b[0], &br, &dzero, &c[0], &m);
  }
  return c;
}

template<> zgrid&
matmult(const zgrid& a, const zgrid& b, zgrid& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw zgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || &c == &b)
    { zgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    zgemm_(&ta, &tb, &m, &n, &k, &zone, &a[0], &ar,
           &b[0], &br, &zzero, &c[0], &m);
  }
  return c;
}

template<> dgrid&
matmult(const dgrid& a, const dvector& b, dgrid& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.size(), bc = 1;
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw dgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || (const dvector*)&c == &b)
    { dgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    dgemm_(&ta, &tb, &m, &n, &k, &done, &a[0],
           &ar, &b[0], &br, &dzero, &c[0], &m);
  }
  return c;
}

template<> zgrid&
matmult(const zgrid& a, const zvector& b, zgrid& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.size(), bc = 1;
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw zgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || (const zvector*)&c == &b)
    { zgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    zgemm_(&ta, &tb, &m, &n, &k, &zone, &a[0],
           &ar, &b[0], &br, &zzero, &c[0], &m);
  }
  return c;
}

template<> dgrid&
matmult(const dvector& a, const dgrid& b, dgrid& c, char ta, char tb)
{
  const int ar = a.size(), ac = 1, br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw dgrid::grid_error("matmult: Grids are not comformable");
  if ((const dvector*)&c == &a || &c == &b)
    { dgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    dgemm_(&ta, &tb, &m, &n, &k, &done, &a[0],
           &ar, &b[0], &br, &dzero, &c[0], &m);
  }
  return c;
}

template<> zgrid&
matmult(const zvector& a, const zgrid& b, zgrid& c, char ta, char tb)
{
  const int ar = a.size(), ac = 1, br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw zgrid::grid_error("matmult: Grids are not comformable");
  if ((const zvector*)&c == &a || &c == &b)
    { zgrid tmp(m, n); c << matmult(a, b, tmp, ta, tb); }
  else {
    c.resize(m, n); c.fill(0);
    zgemm_(&ta, &tb, &m, &n, &k, &zone, &a[0],
           &ar, &b[0], &br, &zzero, &c[0], &m);
  }
  return c;
}

template<> dvector&
matmult(const dgrid& a, const dvector& b, dvector& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.size(), bc = 1;
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw dgrid::grid_error("matmult: Grids are not comformable");
  if (&c == (const dvector*)&a || &c == &b)
    { dvector tmp(m); c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(m); fill(c.begin(), c.end(), 0);
    dgemm_(&ta, &tb, &m, &n, &k, &done, &a[0],
           &ar, &b[0], &br, &dzero, &c[0], &m);
  }
  return c;
}

template<> zvector&
matmult(const zgrid& a, const zvector& b, zvector& c, char ta, char tb)
{
  const int ar = a.rows(), ac = a.cols(), br = b.size(), bc = 1;
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw zgrid::grid_error("matmult: Grids are not comformable");
  if (&c == (const zvector*)&a || &c == &b)
    { zvector tmp(m); c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(m); fill(c.begin(), c.end(), 0);
    zgemm_(&ta, &tb, &m, &n, &k, &zone, &a[0],
           &ar, &b[0], &br, &zzero, &c[0], &m);
  }
  return c;
}

template<> dvector&
matmult(const dvector& a, const dgrid& b, dvector& c, char ta, char tb)
{
  const int ar = a.size(), ac = 1, br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw dgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || &c == (const dvector*)&b)
    { dvector tmp(n); c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(n); fill(c.begin(), c.end(), 0);
    dgemm_(&ta, &tb, &m, &n, &k, &done, &a[0],
           &ar, &b[0], &br, &dzero, &c[0], &m);
  }
  return c;
}

template<> zvector&
matmult(const zvector& a, const zgrid& b, zvector& c, char ta, char tb)
{
  const int ar = a.size(), ac = 1, br = b.rows(), bc = b.cols();
  int m, n, k, k2;
  if (ta == 'n' || ta == 'N') { m = ar; k = ac; } else { k = ar; m = ac; }
  if (tb == 'n' || tb == 'N') { k2 = br; n = bc; } else { n = br; k2 = bc; }

  if (k != k2)
    throw zgrid::grid_error("matmult: Grids are not comformable");
  if (&c == &a || &c == (const zvector*)&b)
    { zvector tmp(n); c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(n); fill(c.begin(), c.end(), 0);
    zgemm_(&ta, &tb, &m, &n, &k, &zone, &a[0],
           &ar, &b[0], &br, &zzero, &c[0], &m);
  }
  return c;
}

//============================================================================
// Compute the eigenvalues, or the magnitudes of the eigenvalues, of a
// square matrix 'a'.
// These routines reasonably handle the case when either the 'evec'
// or 'eval' matrices are the same matrix as 'a' or each other.

template<> int eigenvalues(const dgrid& a, dgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (&eval == &a) {
    dgrid tmp; info = eigenvalues(a, tmp, gridtype, ul); eval << tmp;
    return info;
  }
  if (a.rows() != a.cols()) throw dgrid::grid_error("eigenvalues: Grid is not square");
  if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
  if (gridtype == gengrid) {
    dgrid evec = a; char jobvl = 'N', jobvr = 'N';
    const int n = a.rows(), lwork = 8 * a.rows(); int info = 0;
    eval.resize(n); dgrid evali(n); dvector work(lwork);
    dgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &eval[0], &evali[0],
           0, &n, 0, &n, &work[0], &lwork, &info);
    if (info == 0) {
      for (uint i = 0; i < (uint)n; ++i)
        eval[i] = std::sqrt(eval[i] * eval[i] + evali[i] * evali[i]);
      eval.sort();
    }
  } else if (gridtype == symgrid) {
    dgrid evec = a; char jobz = 'N';
    const int n = a.rows(), lwork = 4 * a.rows();
    eval.resize(n); dvector work(lwork);
    dsyev_(&jobz, &ul, &n, &evec[0], &n, &eval[0], &work[0], &lwork, &info);
  } else
    throw dgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

template<> int eigenvalues(const zgrid& a, dgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (a.rows() != a.cols()) throw zgrid::grid_error("eigenvalues: Grid is not square");
  if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
  if (gridtype == gengrid) {
    zgrid evec = a; char jobvl = 'N', jobvr = 'N';
    const int n = a.rows(), lwork = 4 * a.rows(); int info = 0;
    eval.resize(n); zvector work(lwork); dgrid rwork(2*n);
    zgrid evalz(n);
    zgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &evalz[0],
           0, &n, 0, &n, &work[0], &lwork, &rwork[0], &info);
    for (uint i = 0; i < (uint)n; ++i)
      eval[i] = std::sqrt(std::abs(evalz[i]*evalz[i]));
    eval.sort();
  } else if (gridtype == symgrid) {
    // Need evec, not 'a', since it is destroyed by dheev_
    zgrid evec = a; char jobz = 'N';
    const int n = a.rows(), lwork = 4 * a.rows();
    zvector work(lwork), work2(3 * n - 2);
    eval.resize(n);
    zheev_(&jobz, &ul, &n, &evec[0], &n, &eval[0], &work[0], &lwork,
           &work2[0], &info);
  } else
    throw zgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

template<> int eigenvalues(const dgrid& a, zgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (a.rows() != a.cols()) throw dgrid::grid_error("eigenvalues: Grid is not square");
  if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
  if (gridtype == gengrid) {
    dgrid evec = a; char jobvl = 'N', jobvr = 'N';
    const int n = a.rows(), lwork = 8 * a.rows(); int info = 0;
    eval.resize(n); dgrid evalr(n), evali(n); dvector work(lwork);
    dgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &evalr[0], &evali[0],
           0, &n, 0, &n, &work[0], &lwork, &info);
    if (info == 0)
      for (uint i = 0; i < (uint)n; ++i)
        eval[i] = cplx(evalr[i], evali[i]);
  } else if (gridtype == symgrid) {
    dgrid evec = a; char jobz = 'N';
    const int n = a.rows(), lwork = 4 * a.rows();
    dgrid evald(n); dvector work(lwork);
    dsyev_(&jobz, &ul, &n, &evec[0], &n, &evald[0], &work[0], &lwork, &info);
    if (info == 0) eval = evald;
  } else
    throw dgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

template<> int eigenvalues(const zgrid& a, zgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (&eval == &a) {
    zgrid tmp; info = eigenvalues(a, tmp, gridtype, ul); eval << tmp;
    return info;
  }
  if (a.rows() != a.cols()) throw zgrid::grid_error("eigenvalues: Grid is not square");
  if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
  if (gridtype == gengrid) {
    zgrid evec = a; char jobvl = 'N', jobvr = 'N';
    const int n = a.rows(), lwork = 4 * a.rows(); int info = 0;
    eval.resize(n); zvector work(lwork); dgrid rwork(2*n);
    zgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &eval[0],
           0, &n, 0, &n, &work[0], &lwork, &rwork[0], &info);
  } else if (gridtype == symgrid) {
    // Need 'evec', not 'a', since it is destroyed by dheev_
    zgrid evec = a; char jobz = 'N';
    const int n = a.rows(), lwork = 4 * a.rows();
    zvector work(lwork), work2(3 * n - 2); dgrid evald(n);
    zheev_(&jobz, &ul, &n, &evec[0], &n, &evald[0], &work[0], &lwork,
           &work2[0], &info);
    eval = evald;
  } else
    throw zgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

// Compute the eigenvalues and eigenvectors of a symmetric or hermetian matrix.

template<> int
eigenvectors(const dgrid& a, dgrid& evec, dgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (&eval == &evec)
    throw dgrid::grid_error("eigenvectors: evec and eval are the same matrix");
  if (&eval == &a) {
    dgrid tmp;
    info = eigenvectors(a, evec, tmp, gridtype, ul); eval << tmp;
  } else if (gridtype == gengrid) {
    throw dgrid::grid_error("eigenvectors: Not implemented for general grids");
    dgrid evec = a; char jobvl = 'N', jobvr = 'V';
    const int n = a.rows(), lwork = 8 * a.rows(); int info = 0;
    eval.resize(n); dgrid evali = eval;
    dvector work(lwork);
    dgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &eval[0], &evali[0],
           0, &n, 0, &n, &work[0], &lwork, &info);
  } else if (gridtype == symgrid) {
    evec = a;
    if (a.rows() != a.cols()) throw dgrid::grid_error("eigenvectors: Grid is not square");
    if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
    char jobz = 'V';
    const int n = a.rows(), lwork = 4 * a.rows();
    eval.resize(n); dvector work(lwork);
    dsyev_(&jobz, &ul, &n, &evec[0], &n, &eval[0], &work[0], &lwork, &info);
  } else
    throw dgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

template<> int
eigenvectors(const zgrid& a, zgrid& evec, dgrid& eval, int gridtype, char ul)
{
  int info = 0;
  if (gridtype == gengrid) {
    throw zgrid::grid_error("eigenvectors: Not implemented for general grids");
    zgrid evec = a; char jobvl = 'N', jobvr = 'V';
    const int n = a.rows(), lwork = 4 * a.rows(); int info = 0;
    eval.resize(n); zvector work(lwork); dgrid rwork(2*n);
    zgrid evalz(n);
    zgeev_(&jobvl, &jobvr, &n, &evec[0], &n, &evalz[0],
           0, &n, 0, &n, &work[0], &lwork, &rwork[0], &info);
    for (uint i = 0; i < (uint)n; ++i)
      eval[i] = std::sqrt(std::abs(evalz[i]*evalz[i]));
  } else if (gridtype == symgrid) {
    evec = a;
    if (a.rows() != a.cols()) throw zgrid::grid_error("eigenvectors: Grid is not square");
    if (a.rows() == 0 || a.cols() == 0) { eval.clear(); return(0); }
    char jobz = 'V';
    const int n = a.rows(), lwork = 4 * a.rows();
    zvector work(lwork), work2(3 * n - 2);
    eval.resize(n);
    zheev_(&jobz, &ul, &n, &evec[0], &n, &eval[0], &work[0], &lwork,
           &work2[0], &info);
  } else
    throw zgrid::grid_error("eigenvalues: Unknown grid type");
  return info;
}

//============================================================================
// Compute the singular values of a symmetric or hermetian matrix
// These routines reasonably handle the case when either the 'evec'
// or 'eval' matrices are the same matrix as 'a' or each other.

template<> int singularvalues(const dgrid& a, dgrid& sval)
{
  int info = 0;
  throw dgrid::grid_error("singularvalues: Function not yet implemented");
  return info;
}

template<> int singularvalues(const zgrid& a, dgrid& sval)
{
  int info = 0;
  throw zgrid::grid_error("singularvalues: Function not yet implemented");
  return info;
}

template<> int singularvectors(const dgrid& a, dgrid& svec, dgrid& sval)
{
  int info = 0;
  throw dgrid::grid_error("singularvalues: Function not yet implemented");
  return info;
}

template<> int
singularvectors(const zgrid& a, zgrid& svec, dgrid& sval)
{
  int info = 0;
  throw zgrid::grid_error("singularvalues: Function not yet implemented");
  return info;
}

//============================================================================
// Use Lapack routines to compute the LU decomposition of a
// grid using a LU decomposition appropriate for the type of grid.

// For general and symmetric grids there will be be pivots.
// If the grid type is positive definite, then the LU decomposition
// is the Cholesky decomposition, and the parameter 'ul' is used
// to indicate whether the 'U'pper (default) or 'L'ower Cholesky
// factor is computed.  For triangular matrices, the LU decomposition
// is the original matrix.

template<> int LU(const dgrid& a, dgrid& g, ivector& pivots,
                  int gridtype, char ul)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("LU: Grid is not square");
  const int n = a.rows(); int info = 0;
  if (gridtype == gengrid) {
    g = a; pivots.resize(n);
    dgetrf_(&n, &n, &g[0], &n, &pivots[0], &info);
  } else if (gridtype == symgrid) {
    g = a; pivots.resize(n); dgrid work(n);
    dsytrf_(&ul, &n, &g[0], &n, &pivots[0], &work[0], &n, &info);
  } else if (gridtype == posgrid)
    info = chol(a, g, ul);
  else if (gridtype == trigrid)
    g = a;
  else
    throw dgrid::grid_error("LU: Unknown grid type");
  return info;
}

template<> int LU(const zgrid& a, zgrid& g, ivector& pivots,
                  int gridtype, char ul)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("LU: Grid is not square");
  const int n = a.rows(); int info = 0;
  if (gridtype == gengrid) {
    g = a; pivots.resize(n);
    zgetrf_(&n, &n, &g[0], &n, &pivots[0], &info);
  } else if (gridtype == symgrid) {
    g = a; pivots.resize(n); zgrid work(n);
    zsytrf_(&ul, &n, &g[0], &n, &pivots[0], &work[0], &n, &info);
  } else if (gridtype == posgrid)
    info = chol(a, g, ul);
  else if (gridtype == trigrid)
    g = a;
  else
    throw zgrid::grid_error("LU: Unknown grid type");
  return info;
}

// For positive definite grids, compute the Cholesky Factorization.

template<> int chol(const dgrid& a, dgrid& g, char ul)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("chol: Grid is not square");
  g = a; const int n = (int) a.rows(); int info = 0;
  if (ul == 'U')
    for (int i = 0; i < n - 1; ++i)
      memset(&g[0] + i * n + i + 1, 0, sizeof(double) * (n - i - 1));
  else
    for (int i = 1; i < n; ++i)
      memset(&g[0] + i * n, 0, sizeof(double) * i);
  dpotrf_(&ul, &n, &g[0], &n, &info);
  return info;
}

template<> int chol(const zgrid& a, zgrid& g, char ul)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("chol: Grid is not square");
  g = a; const int n = (int) a.rows(); int info = 0;
  if (ul == 'U')
    for (int i = 0; i < n - 1; ++i)
      memset(&g[0] + i * n + i + 1, 0, sizeof(cplx) * (n - i - 1));
  else
    for (int i = 1; i < n; ++i)
      memset(&g[0] + i * n, 0, sizeof(cplx) * i);
  zpotrf_(&ul, &n, &g[0], &n, &info);
  return info;
}

//============================================================================
// Solve the linear system a * x = b for square 'a' and arbitrary 'b'
// assuming that 'a' contains a previously computed LU decomposition
// of 'a'.  The pivots matrix is ignored in the cases it is not needed.

// Works fine if 'a' or 'b' are the same matrix as 'sol', in which case
// 'a' or 'b' get overwritten with the solution.

template<> int LUsolve(const dgrid& a, const dgrid& b, dgrid& sol,
  const ivector& pivots, int gridtype, char ul, char tr, char dg)
{
  if (&sol == &a) {
    int info = 0; dgrid tmp = b;
    info = LUsolve(a, b, tmp, pivots, gridtype, ul, tr, dg); sol << tmp;
    return info;
  }
  if (a.rows() != a.cols()) throw dgrid::grid_error("LUsolve: Grid is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("LUsolve: Grids are not comformable");
  sol = b; const int n = a.rows(), bc = b.cols(); int info = 0;
  if (gridtype == gengrid)
    dgetrs_(&tr, &n, &bc, &a[0], &n, &pivots[0], &sol[0], &n, &info);
  else if (gridtype == symgrid)
    dsytrs_(&ul, &n, &bc, &a[0], &n, &pivots[0], &sol[0], &n, &info);
  else if (gridtype == posgrid)
    dpotrs_(&ul, &n, &bc, &a[0], &n, &sol[0], &n, &info);
  else if (gridtype == trigrid)
    dtrtrs_(&ul, &tr, &dg, &n, &bc, &a[0], &n, &sol[0], &n, &info);
  else
    throw dgrid::grid_error("LUsolve: Unknown grid type");
  return info;
}

template<> int LUsolve(const zgrid& a, const zgrid& b, zgrid& sol,
  const ivector& pivots, int gridtype, char ul, char tr, char dg)
{
  if (&sol == &a) {
    int info = 0; zgrid tmp = b;
    info = LUsolve(a, b, tmp, pivots, gridtype, ul, tr, dg); sol << tmp;
    return info;
  }
  if (a.rows() != a.cols()) throw zgrid::grid_error("LUsolve: Grid is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("LUsolve: Grids are not comformable");
  sol = b; const int n = a.rows(), bc = b.cols(); int info = 0;
  if (gridtype == gengrid)
    zgetrs_(&tr, &n, &bc, &a[0], &n, &pivots[0], &sol[0], &n, &info);
  else if (gridtype == symgrid)
    zsytrs_(&ul, &n, &bc, &a[0], &n, &pivots[0], &sol[0], &n, &info);
  else if (gridtype == posgrid)
    zpotrs_(&ul, &n, &bc, &a[0], &n, &sol[0], &n, &info);
  else if (gridtype == trigrid)
    ztrtrs_(&ul, &tr, &dg, &n, &bc, &a[0], &n, &sol[0], &n, &info);
  else
    throw zgrid::grid_error("LUsolve: Unknown grid type");
  return info;
}

template<> int cholupdate(dgrid& a, dvector& x)
{
  double ci, si, temp;
  const int n = a.rows();
  for (int i = 0; i < n - 1; ++i) {
    dlartg_(&a[i * (1 + n)], &x[i], &ci, &si, &temp);
    a[i * (1 + n)] = temp;
    int j = n - i - 1;
    drot_(&j, &a[i + (i + 1) * n], &n, &x[i + 1], &ione, &ci, &si);
  }
  dlartg_(&a[n * n - 1], &x[n - 1], &ci, &si, &temp);
  a[n * n - 1] = temp;
  return 0;
}

template<> int cholupdate(zgrid& a, zvector& x)
{
  double ci;
  cplx si, temp;
  const int n = a.rows();
  for (int i = 0; i < n - 1; ++i) {
    zlartg_(&a[i * (1 + n)], &x[i], &ci, &si, &temp);
    a[i * (1 + n)] = temp;
    int j = n - i - 1;
    zrot_(&j, &a[i + (i + 1) * n], &n, &x[i + 1], &ione, &ci, &si);
  }
  zlartg_(&a[n * n - 1], &x[n - 1], &ci, &si, &temp);
  a[n * n - 1] = temp;
  return 0;
}

template<> int cholsolve(const dgrid& a, const dgrid& b, dgrid& sol, char ul)
{
  if (&sol == &a) {
    int info = 0; dgrid tmp = b;
    info = cholsolve(a, b, tmp, ul); sol << tmp; return info;
  }
  if (a.rows() != a.cols()) throw dgrid::grid_error("cholsolve: Grid is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("cholsolve: Grids are not comformable");
  sol = b; const int ar = a.rows(), bc = b.cols(); int info = 0;
  dpotrs_(&ul, &ar, &bc, &a[0], &ar, &sol[0], &ar, &info);
  return info;
}

template<> int cholsolve(const zgrid& a, const zgrid& b, zgrid& sol, char ul)
{
  if (&sol == &a) {
    int info = 0; zgrid tmp = b;
    info = cholsolve(a, b, tmp, ul); sol << tmp; return info;
  }
  if (a.rows() != a.cols()) throw zgrid::grid_error("cholsolve: Grid is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("cholsolve: Grids are not comformable");
  sol = b; const int ar = a.rows(), bc = b.cols(); int info = 0;
  zpotrs_(&ul, &ar, &bc, &a[0], &ar, &sol[0], &ar, &info);
  return info;
}

template<> int trsolve(const dgrid& a, const dgrid& b, dgrid& sol,
                              char ul, char tr, char dg)
{
  if (&sol == &a) {
    int info = 0; dgrid tmp = b;
    info = trsolve(a, b, tmp, ul, tr, dg); sol << tmp; return info;
  }
  if (a.rows() != a.cols()) throw dgrid::grid_error("trsolve: Grid is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("trsolve: Grids are not comformable");
  sol = b; const int ar = a.rows(), bc = b.cols(); int info = 0;
  dtrtrs_(&ul, &tr, &dg, &ar, &bc, &a[0], &ar, &sol[0], &ar, &info);
  return info;
}

template<> int trsolve(const zgrid& a, const zgrid& b, zgrid& sol,
                              char ul, char tr, char dg)
{
  if (&sol == &a) {
    int info = 0; zgrid tmp = b;
    info = trsolve(a, b, tmp, ul, tr, dg); sol << tmp; return info;
  }
  if (a.rows() != a.cols()) throw zgrid::grid_error("trsolve: Grid is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("trsolve: Grids are not comformable");
  sol = b; const int ar = a.rows(), bc = b.cols(); int info = 0;
  ztrtrs_(&ul, &tr, &dg, &ar, &bc, &a[0], &ar, &sol[0], &ar, &info);
  return info;
}

//============================================================================
// Use Lapack routines to compute the inverse of a grid using a LU
// factorization appropriate for the type of grid.

// The paramter 'ul' is used for symmetric and triangular grids and
// indicates whether the 'U'pper (default) or 'L'ower triangle of the
// grid should be used in the computations.  For these matrices, if the
// parameter 'fill' is true (default), then the entire inverse will be
// computed.  Otherwise, only the corresponding upper or lower half will
// be computed, and the other half will essentially be garbage.
// For triangular matrices, 'dg' can be 'N' (default) to indicate the
// diagonal entries will be used, or 'U' to indicate they will be
// assumed to be 1 and ignored.

// Syntax 1: inv(a, ainv)
template<> int inv(const dgrid& a, dgrid& g, int gridtype,
                   char ul, bool fill, char dg)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("inv: Grid is not square");
  g = a; int n = g.rows(); int info = 0;
  if (gridtype == gengrid) {
    ivector ipiv(n); dvector work(n);
    dgetrf_(&n, &n, &g[0], &n, &ipiv[0], &info);
    if (info == 0)
      dgetri_(&n, &g[0], &n, &ipiv[0], &work[0], &n, &info);
  } else if (gridtype == symgrid) {
    ivector ipiv(n); dvector work(n);
    dsytrf_(&ul, &n, &g[0], &n, &ipiv[0], &work[0], &n, &info);
    if (info == 0)
      dsytri_(&ul, &n, &g[0], &n, &ipiv[0], &work[0], &info);
  } else if (gridtype == posgrid) {
    dpotrf_(&ul, &n, &g[0], &n, &info);
    if (info == 0)
      dpotri_(&ul, &n, &g[0], &n, &info);
  } else if (gridtype == trigrid) {
    dtrtri_(&ul, &dg, &n, &g[0], &n, &info);
  }
  if (info == 0 && fill && (gridtype == posgrid || gridtype == symgrid))
    if (ul == 'U')
      for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j) g(i, j) = g(j, i);
    else
      for (int i = 0; i < n - 1; ++i)
       for (int j = i + 1; j < n; ++j) g(i, j) = g(j, i);
  return info;
}

// Syntax 2: ainv = inv(a)
template<> dgrid inv(const dgrid& a, int gridtype, char ul, bool fill, char dg)
  { dgrid tmp; inv(a, tmp, gridtype, ul, fill, dg); return tmp; }

// Syntax 1: inv(a, ainv)
template<> int inv(const zgrid& a, zgrid& g, int gridtype, char ul,
                   bool fill, char dg)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("inv: Grid is not square");
  g = a; const int n = g.rows(); int info = 0;
  if (gridtype == gengrid) {
    ivector ipiv(n); zvector work(n);
    zgetrf_(&n, &n, &g[0], &n, &ipiv[0], &info);
    if (info == 0)
      zgetri_(&n, &g[0], &n, &ipiv[0], &work[0], &n, &info);
  } else if (gridtype == posgrid) {
    zpotrf_(&ul, &n, &g[0], &n, &info);
    if (info == 0)
      zpotri_(&ul, &n, &g[0], &n, &info);
  } else if (gridtype == symgrid) {
    ivector ipiv(n); zvector work(n);
    zsytrf_(&ul, &n, &g[0], &n, &ipiv[0], &work[0], &n, &info);
    if (info == 0)
      zsytri_(&ul, &n, &g[0], &n, &ipiv[0], &work[0], &info);
  } else if (gridtype == trigrid) {
    ztrtri_(&ul, &dg, &n, &g[0], &n, &info);
  }
  if (info == 0 && fill && (gridtype == posgrid || gridtype == symgrid))
    if (ul == 'U')
      for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j) g(i, j) = g(j, i);
    else
      for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j) g(i, j) = g(j, i);
  return info;
}

// Syntax 2: ainv = a.inv()
template<> zgrid inv(const zgrid& a, int gridtype, char ul,
                     bool fill, char dg)
  { zgrid tmp; inv(a, tmp, gridtype, ul, fill, dg); return tmp; }

//============================================================================
// Compute y = alpha * x + y

dgrid& axpy(const dgrid& x, double alpha, dgrid& y)
{
  if (x.size() != y.size()) throw dgrid::grid_error("axpy: Grids are not comformable");
  const int as = x.size();
  if (as != 0) daxpy_(&as, &alpha, &x[0], &ione, &y[0], &ione);
  return y;
}

zgrid& axpy(const zgrid& x, cplx alpha, zgrid& y)
{
  if (x.size() != y.size()) throw zgrid::grid_error("axpy: Grids are not comformable");
  const int as = x.size();
  if (as != 0) zaxpy_(&as, &alpha, &x[0], &ione, &y[0], &ione);
  return y;
}

//============================================================================
// Compute c = alpha * a * b + beta * c

dgrid& trmm(const dgrid& a, dgrid& b, double alpha, const char side,
            const char uplo, const char tra, const char diag)
{
  if (a.cols() != b.rows())
    throw dgrid::grid_error("trmm: Grids are not comformable");
  if (&a == &b)
    throw dgrid::grid_error("trmm: The two grid arguments are the same grid");
  if (a.size() != 0 && b.size() != 0) {
    const int ar = a.rows(), ac = a.cols(), bc = b.cols();
    dtrmm_(&side, &uplo, &tra, &diag, &ac, &bc, &alpha, &a[0], &ar, &b[0], &ac);
  }
  return b;
}

zgrid& trmm(const zgrid& a, zgrid& b, const cplx& alpha, const char side,
            const char uplo, const char tra, const char diag)
{
  if (a.cols() != b.rows())
    throw zgrid::grid_error("trmm: Grids are not comformable");
  if (&a == &b)
    throw zgrid::grid_error("trmm: The two grid arguments are the same grid");
  if (a.size() != 0 && b.size() != 0) {
    const int ar = a.rows(), ac = a.cols(), bc = b.cols();
    ztrmm_(&side, &uplo, &tra, &diag, &ac, &bc, &alpha, &a[0], &ar, &b[0], &ac);
  }
  return b;
}

dgrid& gemm(const dgrid& a, const dgrid& b, dgrid& c, double alpha,
            double beta, char tra, char trb)
{
  if (a.cols() != b.rows() || a.rows() != c.rows() || b.cols() != c.cols())
    throw dgrid::grid_error("gemm: Grids are not comformable");
  if (&a == &c || &b == &c)
    throw dgrid::grid_error("gemm: Result grid is the same as one of the multipliers");
  if (a.cols() != 0 && c.size() != 0) {
    const int ar = a.rows(), ac = a.cols(), bc = b.cols();
    dgemm_(&tra, &trb, &ar, &bc, &ac, &alpha, &a[0], &ar,
           &b[0], &ac, &beta, &c[0], &ar);
  }
  return c;
}

zgrid& gemm(const zgrid& a, const zgrid& b, zgrid& c, const cplx& alpha,
            const cplx& beta, char tra, char trb)
{
  if (a.cols() != b.rows() || a.rows() != c.rows() || b.cols() != c.cols())
    throw zgrid::grid_error("gemm: Grids are not comformable");
  if (&a == &c || &b == &c)
    throw zgrid::grid_error("gemm: Result grid is the same as one of the multipliers");
  if (a.cols() != 0 && c.size() != 0) {
    const int ar = a.rows(), ac = a.cols(), bc = b.cols();
    zgemm_(&tra, &trb, &ar, &bc, &ac, &alpha, &a[0], &ar,
           &b[0], &ac, &beta, &c[0], &ar);
  }
  return c;
}

//============================================================================
// Solve the over/under-determined linear system a * x = b for arbitary
// 'a' and 'b' using a QR or LQ factorization.  On exit, 'a' contains
// its QR or LQ factorization, while the columns of 'b' contain the
// solution 'x'.

int gels(dgrid& a, dgrid& b, char trans)
{
  if (a.rows() != b.rows()) throw dgrid::grid_error("gels: Grids are not comformable");
  const int m = a.rows(), n = a.cols(), nrhs = b.cols(); int info = 0;
  const int lwork = 10 * (m + n + nrhs);
  dvector work(lwork);
  dgels_(&trans, &m, &n, &nrhs, &a[0], &m, &b[0], &m,
         &work[0], &lwork, &info);
  return info;
}

int gels(zgrid& a, zgrid& b, char trans)
{
  if (a.rows() != b.rows()) throw zgrid::grid_error("gels: Grids are not comformable");
  const int m = a.rows(), n = a.cols(), nrhs = b.cols(); int info = 0;
  const int lwork = 10 * (m + n + nrhs);
  zvector work(lwork);
  zgels_(&trans, &m, &n, &nrhs, &a[0], &m, &b[0], &m,
         &work[0], &lwork, &info);
  return info;
}

//============================================================================
// Solve the linear system a * x = b for square 'a' and arbitrary 'b'.
// On exit, if 'a' and 'b' are not the same matrix, 'a' contains an LU or
// Cholesky decomposition, and 'b' contains the solution 'x'.  If 'a' and
// 'b' are the same matrix, then they will contain the identity solution.

// General square 'a'

int gesv(dgrid& a, dgrid& b, ivector& pivots)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("gesv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("gesv: Grids are not comformable");
  if (&a == &b) { b = dgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0; pivots.resize(n);
  dgesv_(&n, &nrhs, &a[0], &n, &pivots[0], &b[0], &n, &info);
  return info;
}

int gesv(zgrid& a, zgrid& b, ivector& pivots)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("gesv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("gesv: Grids are not comformable");
  if (&a == &b) { b = zgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0; pivots.resize(n);
  zgesv_(&n, &nrhs, &a[0], &n, &pivots[0], &b[0], &n, &info);
  return info;
}

// Symmetric square 'a'

int sysv(dgrid& a, dgrid& b, ivector& pivots, char ul)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("sysv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("sysv: Grids are not comformable");
  if (&a == &b) { b = dgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0; pivots.resize(n);
  const int lwork = 2 * n; dgrid work(lwork);
  dsysv_(&ul, &n, &nrhs, &a[0], &n, &pivots[0],
         &b[0], &n, &work[0], &lwork, &info);
  return info;
}

int sysv(zgrid& a, zgrid& b, ivector& pivots, char ul)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("sysv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("sysv: Grids are not comformable");
  if (&a == &b) { b = zgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0; pivots.resize(n);
  const int lwork = 2 * n; zgrid work(lwork);
  zsysv_(&ul, &n, &nrhs, &a[0], &n, &pivots[0],
         &b[0], &n, &work[0], &lwork, &info);
  return info;
}

// Hermetian square 'a'

int hesv(zgrid& a, zgrid& b, ivector& pivots, char ul)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("hesv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("hesv: Grids are not comformable");
  if (&a == &b) { b = zgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0; pivots.resize(n);
  const int lwork = 2 * n; zgrid work(lwork);
  zhesv_(&ul, &n, &nrhs, &a[0], &n, &pivots[0],
         &b[0], &n, &work[0], &lwork, &info);
  return info;
}

// Positive Definite square 'a'

int posv(dgrid& a, dgrid& b, char ul)
{
  if (a.rows() != a.cols()) throw dgrid::grid_error("posv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw dgrid::grid_error("posv: Grids are not comformable");
  if (&a == &b) { b = dgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0;
  dposv_(&ul, &n, &nrhs, &a[0], &n, &b[0], &n, &info);
  return info;
}

int posv(zgrid& a, zgrid& b, char ul)
{
  if (a.rows() != a.cols()) throw zgrid::grid_error("posv: Grid 'A' is not square");
  if (a.rows() != b.rows()) throw zgrid::grid_error("posv: Grids are not comformable");
  if (&a == &b) { b = zgrid("I", a.cols()); return 0; }
  const int n = a.rows(), nrhs = b.cols(); int info = 0;
  zposv_(&ul, &n, &nrhs, &a[0], &n, &b[0], &n, &info);
  return info;
}

//============================================================================
// Solve a triangular system (*this) * x = b for x using either
// back or forward substitution.

template<> dgrid&
backsolve(const dgrid& a, const dgrid& b, dgrid& x, char ul) {
  x = b;
  if (ul == 'U') trsm(a, b, x, 1.0, 'L', 'U', 'N', 'N');
  else trsm(a, b, x, 1.0, 'L', 'L', 'T', 'N');
  return x;
}

template<> zgrid&
backsolve(const zgrid& a, const zgrid& b, zgrid& x, char ul) {
  x = b;
  if (ul == 'U') trsm(a, b, x, 1.0, 'L', 'U', 'N', 'N');
  else trsm(a, b, x, 1.0, 'L', 'L', 'T', 'N');
  return x;
}

template<> dgrid&
forwardsolve(const dgrid& a, const dgrid& b, dgrid& x, char ul) {
  x = b;
  if (ul == 'L') trsm(a, b, x, 1.0, 'L', 'L', 'N', 'N');
  else trsm(a, b, x, 1.0, 'L', 'U', 'T', 'N');
  return x;
}

template<> zgrid&
forwardsolve(const zgrid& a, const zgrid& b, zgrid& x, char ul) {
  x = b;
  if (ul == 'L') trsm(a, b, x, 1.0, 'L', 'L', 'N', 'N');
  else trsm(a, b, x, 1.0, 'L', 'U', 'T', 'N');
  return x;
}

//============================================================================
// Solve prod(op(a), x) = alpha * b, where 'a' is a triangular matrix.
// If 'sd' == 'L', prod(a,b) = a %*% b, otherwise, prod(a,b) = b %*% a.
// If 'ul' == 'U', use upper half of 'a', otherwise use lower half.
// If 'tr' == 'N', op(a) = a, otherwise op(a) = t(a).
// If 'dg' != 'U', assume diagonal entries of 'a' are 1.

// If x and b are equal, b will be overwritten with the solution x.

dgrid& trsm(const dgrid& a, const dgrid& b, dgrid& x, double alpha,
            char sd, char ul, char tr, char dg)
{
  if (sd != 'l' && sd != 'L') sd = 'R';
  if (ul != 'u' && ul != 'U') ul = 'L';
  if (tr != 'n' && tr != 'N') tr = 'T';
  if (dg != 'n' && dg != 'N') dg = 'U';
  x = b; const int br = b.rows(), bc = b.cols(), ar = a.rows();
  dtrsm_(&sd, &ul, &tr, &dg, &br, &bc, &alpha, &a[0], &ar, &x[0], &br);
  return x;
}

zgrid& trsm(const zgrid& a, const zgrid& b, zgrid& x, const cplx& alpha,
            char sd, char ul, char tr, char dg)
{
  if (sd != 'l' && sd != 'L') sd = 'R';
  if (ul != 'u' && ul != 'U') ul = 'L';
  if (tr != 'n' && tr != 'N') tr = 'T';
  if (dg != 'n' && dg != 'N') dg = 'U';
  x = b; const int br = b.rows(), bc = b.cols(), ar = a.rows();
  ztrsm_(&sd, &ul, &tr, &dg, &br, &bc, &alpha, &a[0], &ar, &x[0], &br);
  return x;
}

//============================================================================
// grid.C

#endif