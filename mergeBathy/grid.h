/** 
* @file			grid.h
* @brief		Grid class used throughout mergeBathy.
*
*/

#pragma once
// ==========================================================================
// Issues:

// Should sorting work for complex ??  There is a dummy template that
// doesn't do anything.

//Downloaded from http:\\oldmill.uchicago.edu/~wilder/Code/grid/

// ==========================================================================
// Includes
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <complex>
#include <vector>
#include <iterator>
#include <functional>
#include <numeric>
#include <algorithm>
#include <string>
#include <cstring>

#ifndef USE_ASSERTS
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

using std::vector; using std::string;

#ifndef uchar
typedef unsigned char uchar;
#endif
#ifndef uint
typedef unsigned int uint;
#endif

typedef std::complex<double> cplx;

enum { gengrid = 0, posgrid = 1, symgrid = 2, trigrid = 3};

// ==========================================================================
// Grid (matrix) of elements of type T.
// A grid is essentially a standard library vector with added functionality.

template<class T>
class grid
{
	protected:
		vector<T> sto;
		uint nr, nc; // number of rows and columns

		// Consistency Checking Functions
		bool inrange(int r, int c) const 
			{ return (r >= 0 && c >= 0 && r < (int)nr && c < (int)nc); }
		bool inrange(uint r, uint c) const { return (r < nr && c < nc); }
		bool inrange1(uint r) const { return (r < (nr * nc)); }
		template<class S> bool inrange(S r, S c) const
			{ return (r < (S) nr && c < (S) nc); }
		template<class S> bool inrange1(S r) const
			{ return (r < (S)(nr * nc)); }

	public:
		typedef typename vector<T>::iterator iterator;
		typedef typename vector<T>::const_iterator const_iterator;
		typedef size_t size_type;
		typedef T value_type;
		typedef T& reference;
		typedef const T& const_reference;

		// Constructors, Destructor, operator= and operator==
		grid() : nr(0), nc(0) { }
		explicit grid(uint r, uint c=1) : sto(r * c), nr(r), nc(c) { }
		grid(uint r, uint c, T v) : sto(r * c, v), nr(r), nc(c) { }
		explicit grid(const string& s) : nr(0), nc(0) { read(s); } // File input
		grid(const string& s, uint N); // Special Matrices
		grid(const string&, T, T, T); // Sequence
		grid(const T* tp, uint r, uint c) : sto(r * c), nr(r), nc(c)
			{ std::copy(tp, tp + r * c, begin()); }
		grid(const vector<T>& v) : sto(v), nr((uint)v.size()), nc(1) { }
		grid(const grid<T>& g) : sto(g.sto), nr(g.nr), nc(g.nc) { }
		~grid() { }
		grid<T>& operator=(const grid<T> &g);
		template<class S> void operator=(const grid<S> &g);
		bool operator==(const grid<T> &g) const;
		operator const vector<T>&() const { return sto; }
		// The following non-const conversion can be dangerous.
		// operator vector<T>&() { return sto; }
		vector<T>& vec() { return sto; }

		// rows(), cols()
		uint rows() const { return nr; }
		uint cols() const { return nc; }

		// Basic Member Functions
		void clear() { nr = nc = 0; sto.clear(); }
		uint size() const { return (uint) sto.size(); }
		T& front() { return sto.front(); }
		const T& front() const { return sto.front(); }
		iterator begin() { return sto.begin(); }
		iterator end() { return sto.end(); }
		const_iterator begin() const { return sto.begin(); }
		const_iterator end() const { return sto.end(); }
		void resize(uint r, uint c = 1) { sto.resize(r * c); nr = r; nc = c; }
		void resize(uint r, uint c, const T& v)
			{ sto.resize(r * c, v); nr = r; nc = c; }
		void fill(const T& val = 0) { std::fill(begin(), end(), val); }
		void push_back(const T& t)
			{ if (nc == 1) ++nr; else if (nr == 1) ++nc; else return;
			sto.push_back(t); }

		// Various versions of operator()
		T& operator[](uint r)
			{ assert(inrange1(r)); return sto[r]; }
		template<class S> T& operator[](S r)
			{ assert(inrange1(r)); return sto[r]; }
		template<class S> const T& operator[](S r) const
			{ assert(inrange1(r)); return sto[r]; }

		T& operator()(uint r)
			{ assert(inrange1(r)); return sto[r]; }
		template<class S> T& operator()(S r)
			{ assert(inrange1(r)); return sto[r]; }
		template<class S> const T& operator()(S r) const
			{ assert(inrange1(r)); return sto[r]; }

		T& operator()(uint r, uint c)
			{ assert(inrange(r, c)); return sto[c * nr + r]; }
		const T& operator()(uint r, uint c) const
			{ assert(inrange(r, c)); return sto[c * nr + r]; }
		template<class S> T& operator()(S r, S c)
			{ assert(inrange(r, c)); return sto[c * nr + r]; }
			template<class S> const T& operator()(S r, S c) const
			{ assert(inrange(r, c)); return sto[c * nr + r]; }

		template<class S> T& operator()(const S* p)
			{ assert(inrange(p[0], p[1])); return(operator()(p[0], p[1])); } 
		template<class S> const T& operator()(const S* p) const
			{ assert(inrange(p[0], p[1])); return(operator()(p[0], p[1])); } 

		// Safe operator(), but very inefficient.
		void set(int r, int c, T val)
			{ if (inrange(r, c)) operator()(r, c) = val; }
		T get(int r, int c) const
			{ return (inrange(r, c) ? operator()(r, c) : 0); }

		T max() const 
			{ return (size() == 0) ? 0 : *std::max_element(begin(), end()); }
		T min() const 
			{ return (size() == 0) ? 0 : *std::min_element(begin(), end()); }

		// I/O Functions
		bool write(const string& file, bool header = 1, bool binary = 1);
		bool write(std::ofstream& os, bool binary = 1);
		bool read(const string& file, bool header = 1);
		bool read(std::ifstream& is, bool binary = 1);
		int memory_size() const;

		// Manipulate the contents of grids.
		void swap(vector<T>& v)
			{ nr = v.size(); nc = 1; sto.swap(v); }
		void swap(grid<T>& g) {
			if (this != &g) 
				{ sto.swap(g.sto); std::swap(nr, g.nr); std::swap(nc, g.nc); }
		}
		grid<T>& operator<<(vector<T> &v)
			{ swap(v); v.clear(); return *this; }
		grid<T>& operator<<(grid<T> &g) {
			if (this != &g) { swap(g); g.sto.clear(); g.nr = g.nc = 0; }
				return *this;
		}

		void subgrid(grid<T>& g, uint r, uint c, uint numrows, uint numcols);
		void rbind(const grid<T>& g);
		void cbind(const grid<T>& g);

		// Return the number of 0 or non-0 pixels in a grid.
		uint offpixels() const { 
			uint c = 0;
			for (const_iterator i = begin(); i != end(); ++i)
				c += (*i == 0);
			return c;
		}
		uint onpixels() const { return(size() - offpixels()); }

		// Useful Utility Functions
		void normalize(T val1, T val2);
		grid<T>& sort() { std::sort(begin(), end()); return *this; }

	// Exceptions
	class grid_error : public std::exception {
		private:
			string msg;
		public:
			virtual const char* what() const throw() { return(msg.c_str()); }
			grid_error(const string& s = "") : msg(s) { }
			~grid_error() throw() { }
	};

}; // class grid<T>

// ==========================================================================
// Useful typedefs

typedef vector<bool> bvector;
typedef vector<char> cvector;
typedef vector<unsigned char> ucvector;
typedef vector<int> ivector;
typedef vector<uint> uivector;
typedef vector<long> lvector;
typedef vector<float> fvector;
typedef vector<double> dvector;
typedef vector<cplx> zvector;
typedef vector<string> strvector;

typedef bvector::iterator biter;
typedef cvector::iterator citer;
typedef ucvector::iterator uciter;
typedef ivector::iterator iiter;
typedef uivector::iterator uiiter;
typedef lvector::iterator liter;
typedef fvector::iterator fiter;
typedef dvector::iterator diter;
typedef zvector::iterator ziter;
typedef strvector::iterator striter;

typedef bvector::const_iterator bconstiter;
typedef cvector::const_iterator cconstiter;
typedef ucvector::const_iterator ucconstiter;
typedef ivector::const_iterator iconstiter;
typedef uivector::const_iterator uiconstiter;
typedef lvector::const_iterator lconstiter;
typedef fvector::const_iterator fconstiter;
typedef dvector::const_iterator dconstiter;
typedef zvector::const_iterator zconstiter;
typedef strvector::const_iterator strconstiter;

typedef grid<bool> bgrid;
typedef grid<char> cgrid;
typedef grid<unsigned char> ucgrid;
typedef grid<int> igrid;
typedef grid<uint> uigrid;
typedef grid<long> lgrid;
typedef grid<float> fgrid;
typedef grid<double> dgrid;
typedef grid<cplx> zgrid;
typedef grid<string> strgrid;

typedef vector<bgrid> bgrids;
typedef vector<cgrid> cgrids;
typedef vector<ucgrid> ucgrids;
typedef vector<igrid> igrids;
typedef vector<uigrid> uigrids;
typedef vector<lgrid> lgrids;
typedef vector<fgrid> fgrids;
typedef vector<dgrid> dgrids;
typedef vector<zgrid> zgrids;

typedef vector<bgrids> bgridss;
typedef vector<cgrids> cgridss;
typedef vector<ucgrids> ucgridss;
typedef vector<igrids> igridss;
typedef vector<uigrids> uigridss;
typedef vector<lgrids> lgridss;
typedef vector<fgrids> fgridss;
typedef vector<dgrids> dgridss;
typedef vector<zgrids> zgridss;

//============================================================================

template<> inline grid<cplx>& grid<cplx>::sort() { return *this; }

// ==========================================================================
// Special matrix constructors

template<class T>
grid<T>::grid(const string& s, uint N) : sto(N * N, 0), nr(N), nc(N)
{
  if (s == "I" || s.substr(0, 3) == "ide") {
    // N x N Identity Matrix
    for (uint i = 0; i != N; ++i)
      operator()(i, i) = 1;
  } else if (s.substr(0, 3) == "hil") {
    // N x N Badly conditioned Hilbert Matrix
    for (uint i = 0; i != N; ++i)
      for (uint j = 0; j != N; ++j)
        operator()(i, j) = T(1.0 / (i + j + 1.0));
  } else if (s.substr(0, 3) == "wil") {
    // N x N Wilkinson Eigenvalue Test Matrix
    for (uint i = 0; i != N - 1; ++i) {
      operator()(i, i + 1) = 1;
      operator()(i + 1, i) = 1;
    }
    for (uint i = 0; i != N; ++i)
      operator()(i, i) = T(std::abs(i - (N - 1) / 2.0));
  } else
    throw grid_error("grid(string&, uint): Unknown initialization");
}

template<class T>
grid<T>::grid(const string& s, T start, T end, T step)
{
  if (s.substr(0, 3) == "seq") {
    uint count = (uint) std::abs((end - start) / step) + 1;
    resize(count);
    for (uint i = 0; i != count; ++i) {
      sto[i] = start;
      start += step;
    }
  } else
    throw grid_error("grid(string&, T, T, T): Unknown initialization");
}

// ==========================================================================
// Print a grid; Will only work for classes with ostream<<(const T&).
// If 'tr' is non-zero, the output will be transposed, and if 'max'
// is non-zero, at most 'max' rows will be output.  This function
// provides somewhat finer control over the output than operator<<.
// By default, fields are separated by spaces.  An alternative separator
// can be specified, but this is designed for compatibility with other
// programs as then the file can not be easily read back using this code.

template<class T>
void dump(const grid<T>& g, std::ostream& os=std::cout,
          uint tr=0, uint max=0, int prec=4, char sep=' ')
{
  if (g.size() == 0) { os << "<Empty grid>\n"; return; }
  int wid = prec + 6;
  if (tr) {
    if (max == 0 || max > g.cols()) max = g.cols();
    for (uint j = 0; j != max; ++j) {
      for (uint i = 0; i != g.rows(); ++i)
        os << std::setw(wid) << std::setprecision(prec)
                  << (g(i, j)) << sep;
      os << "\n";
    }
  } else {
    if (max == 0 || max > g.rows()) max = g.rows();
    for (uint i = 0; i != max; ++i) {
      for (uint j = 0; j != g.cols(); ++j)
        os << std::setw(wid) << std::setprecision(prec)
                  << (g(i, j)) << sep;
      os << "\n";
    }
  }
}

// Call the above grid<T>::dump with a fixed field output.

template<class T>
void dumpfixed(const grid<T>& g, std::ostream& os=std::cout,
               uint tr=0, uint max=0, int prec=4, char sep=' ')
{
  
  std::ios::fmtflags f = os.setf(std::ios::fixed);
  dump(g, os, tr, max, prec, sep);
  os.flags(f);
}

// Print a vector/grid to an ostream using a default format similar to
// the one used in matlab.  A newline is printed before the vector/grid.

template<class T>
std::ostream& operator<<(std::ostream& os, const vector<T>& v)
{
  if (v.size() == 0) { os << "<Empty vector>\n"; return os; }
  const int wid = 10, prec = 4;
  std::ios::fmtflags f = os.setf(std::ios::fixed);
  for (uint i = 0; i != v.size(); ++i)
    os << std::setw(wid) << std::setprecision(prec) << v[i] << " ";
  os << "\n";
  os.flags(f);
  return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const grid<T>& g)
  { dumpfixed(g, os, 0, 0, 4); return os; }

// ==========================================================================
// operator<< specialized to unsigned char, defined in grid.C

template<> std::ostream& operator<<(std::ostream& os, const ucgrid& g);

// ==========================================================================
// (Binary) write out a grid to an open output stream.

template<class T>
bool grid<T>::write(std::ofstream& ofs, bool binary)
{
  if (!ofs)
    return false;

  if (binary) {
    ofs.write((char*) &nr, sizeof(nr));
    ofs.write((char*) &nc, sizeof(nc));
    for (uint i = 0; i != nr * nc; ++i)
      ofs.write((char*) &sto[i], sizeof(T));
  } else {
    ofs << nr << " " << nc << "\n";
    for (uint i = 0; i != rows(); ++i) {
      for (uint j = 0; j != cols(); ++j)
        ofs << grid<T>::operator()(i, j);
      ofs << "\n";
    }
  }
  return true;

} // grid::write

// ==========================================================================
// (Binary) write out a grid to a file.

template<class T>
bool grid<T>::write(const string& file, bool header, bool binary)
{
  std::ofstream ofs(file.c_str());
  if (!ofs)
    return false;

  if (header)
    ofs.write("GR11", 4);
  return write(ofs, binary);

} // grid::write


// ==========================================================================
// Calculate the total amount of disk space used by a binary written grid.

template<class T> int grid<T>::memory_size() const
  { return (sizeof(nr) + sizeof(nc) + size() * sizeof(T)); }

// ==========================================================================
// Read in a grid which was previously written to a file

template<class T> bool grid<T>::read(const string& file, bool header)
{
  std::ifstream ifs(file.c_str());
  if (!ifs)
    return false;

  char version[4];
  if (header == 0)
    memcpy(version, "GR12", 4);
  else
    ifs.read(version, 4);
  if (memcmp(version, "GR11", 4) != 0 && memcmp(version, "GR12", 4) != 0) {
    if (loadpgm(*this, file))
      return true;
    else {
      std::cerr << "[" << file << "] is not a grid or pgm file" << std::endl;
      return false;
    }
  }

  // GR11 indicates a platform dependent binary file
  if (!memcmp(version, "GR11", 4))
    read(ifs, 1);
  // GR12 should be a platform independent text file
  else if (!memcmp(version, "GR12", 4))
    read(ifs, 0);
  return true;

} // grid::read

// ==========================================================================
// Read in a grid from an open istream

template<class T>
bool grid<T>::read(std::ifstream& ifs, bool binary)
{
  if (!ifs || ifs.eof()) return false;

  if (binary) {
    ifs.read((char*) &nr, sizeof(nr));
    ifs.read((char*) &nc, sizeof(nc));
    resize(nr, nc);
    for (uint i = 0; i != nr * nc; ++i)
      ifs.read((char*) &sto[i], sizeof(T));
  } else {
    ifs >> nr;
    ifs >> nc;
    resize(nr, nc);
    for (uint i = 0; i != rows(); ++i)
      for (uint j = 0; j != cols(); ++j)
        ifs >> grid<T>::operator()(i, j);
  }
  return true;

} // grid::read

//============================================================================

template<class T>
bool loadpgm(grid<T>& a, const string& pgmname)
{
  // Open the pgm file
  std::ifstream ifs(pgmname.c_str());
  if (!ifs)
    return false;

  // Find the first line that doesn't begin with white space.
  char buf[255];
  char pchar;
  int mode, r, c, maxval, matches;
  do {
    ifs.getline(buf, 255);
  } while (!ifs.eof() && (buf[0] == 0 || buf[0] == ' '));
  if (ifs.eof())
    return false;

  // Make sure the file is a pgm file.  Determine the pgm mode, the
  // dimensions, and the range of pixel values
  matches = sscanf_s(buf, "%c%1d%d%d%d", &pchar, &mode, &r, &c, &matches);
  if (matches < 2 || mode < 2 || mode > 6)
    return false;

  if (matches < 5) {
    ifs.getline(buf, 255);
    while (!ifs.eof() && (buf[0] == '#' || buf[0] == 0))
    ifs.getline(buf, 255);
    sscanf_s(buf, "%d %d", &r, &c);
    ifs.getline(buf, 255);
    sscanf_s(buf, "%d", &maxval);
    if (ifs.eof() || r <= 0 || c <= 0 || maxval <= 0)
      return false;
  }

  a.resize(r, c);
  int i, j;
  switch(mode)
  {
   case 2:
    for(j = 0; j < c; ++j) {
      int tmpi;
      for(i = 0; i < r; ++i) {
        ifs >> tmpi;
        a(i, j) = tmpi;
      }
    };
    break;

   case 3:
    for(j = 0; j < c; j++) {
      int tmpi[3];
      for(i = 0; i < r; ++i) {
        ifs >> tmpi[0] >> tmpi[1] >> tmpi[2];
        a(i, j) = (unsigned char)
          (0.212671 * tmpi[0] + 0.715160 * tmpi[1] + 0.072169 * tmpi[2]);
      }
    };
    break;

   case 5:
    for(j = 0; j < c; ++j) {
      for(i = 0; i < r; ++i)
        ifs.read((char *)&a(i, j), 1);
    };
    break;

   case 6:
    for(j=0; j < c; j++) {
      unsigned char tmpc[3];
      for(i = 0; i < r; ++i) {
        ifs.read((char *)tmpc, 3);
        a(i, j) = (unsigned char)
          (0.212671 * tmpc[0] + 0.715160 * tmpc[1] + 0.072169 * tmpc[2]);
      }
    };
    break;

   default:
    return false;
    break;
  };
  return true;

} // loadpgm

// ==========================================================================
// Set a grid equal to a grid 'g'; leave 'g' unchanged.
// Note: this routine is very efficient.

template<class T> grid<T>& grid<T>::operator=(const grid<T>& g) {
  if (this != &g)
    { resize(g.rows(), g.cols()); std::copy(g.begin(), g.end(), begin()); }
  return *this;
}

template<class T> template<class S>
void grid<T>::operator=(const grid<S>& g) {
  if ((void*)this != (void*)&g)
    { resize(g.rows(), g.cols()); std::copy(g.begin(), g.end(), begin()); }
}

// ==========================================================================
// Test if two matrices are equal.

template<class T>
bool grid<T>::operator==(const grid<T> &g) const
{
  if (this == &g) return true;
  if (size() != g.size()) return false;
  return std::equal(g.begin(), g.end(), begin());
}

// ==========================================================================
// Set 'g' equal to a subgrid of *this.  ('r', 'c') is the starting point
// in 'g', and the dimension of the subgrid is ('rsize', 'csize').

template<class T>
void grid<T>::subgrid(grid<T>& g, uint r, uint c, uint rsize, uint csize)
{
  if (rsize == 0 || csize == 0 || r >= rows() || c >= cols())
    { g.clear(); return; }

  if (rsize > rows() - r) rsize = rows() - r;
  if (csize > cols() - c) csize = cols() - c;
  if (this == &g)
    { grid tmp; subgrid(tmp, r, c, rsize, csize); swap(tmp); }
  else {
    g.resize(rsize, csize);
    T *gp = &g[0], *thisp = &(operator()(r, c));
    for (uint j = 0; j != csize; ++j) {
      for (uint i = 0; i != rsize; ++i)
        *gp++ = *thisp++;
      thisp += rows() - rsize;
    }
  }

} // grid::subgrid

// ==========================================================================
// Append the rows or columns of 'g' to '*this'.  Correctly allows
// 'g' to be appended to itself.

template<class T> void grid<T>::cbind(const grid<T>& g)
{
  if (rows() != g.rows())
    throw grid_error("cbind: Grids are not conformable");
  uint nx = rows(), ny = cols(), ng = g.size();
  resize(nx, ny + g.cols());
  std::copy(g.begin(), g.begin() + ng, begin() + nx * ny);
}

template<class T> void grid<T>::rbind(const grid<T>& g)
{
  if (cols() != g.cols())
    throw grid_error("rbind: Grids are not conformable");
  uint nx = rows(), ny = cols(), gx = g.rows();
  grid<T> tmp(nx + gx, ny);
  T* thisp = &front(), *tmpp = &tmp[0];
  const T* gp = &g(0, 0);
  for (uint i = 0; i != ny; ++i) {
    std::copy(thisp, thisp + nx, tmpp); 
    std::copy(gp, gp + gx, tmpp + nx);
    thisp += nx; gp += gx; tmpp += nx + gx;
  }
  *this << tmp;
}

// ==========================================================================
// Add/Subtract/Multiply/Divide a vector 'v' to/by a vector or grid 'g'.
// Column order, recycle vector elements, so arguments do not need to
// have the same size.

#define VGOP1(VG,OP1,OP2) \
template<class T> VG <T>& \
operator OP1(VG <T>& g, const vector<T>& v) { \
  if (g.size() == 0 || v.size() == 0) return g; \
  T* gp = &g[0]; const T* vp = &v[0]; \
  for (uint i = 0; i != g.size(); ++i) { \
    if (i % v.size() == 0) vp = &v[0]; \
    *gp++ OP1 *vp++; \
  } \
  return g; \
} \
template<class T> VG <T> \
operator OP2(const VG <T>& g, const vector<T>& v) \
  { VG <T> tmp = g; return operator OP1(tmp, v); }

#define VGOP(OP1, OP2) VGOP1(grid, OP1, OP2) VGOP1(vector, OP1, OP2)

VGOP(+=, +) VGOP(-=, -) VGOP(*=, *) VGOP(/=, /)

#undef VGOP1
#undef VGOP

// Add/Subtract/Multiply/Divide a grid 'g2' to/by a grid 'g1'
// For these operations the grids must be conformable.

#define GGOP(OP1,OP2) \
template<class T> grid<T>& \
operator OP1(grid<T>& g1, const grid<T>& g2) { \
  if (g1.rows() != g2.rows() || g1.cols() != g2.cols()) \
    throw typename grid<T>::grid_error("operator" # OP1 ": Grids are not conformable"); \
  T *g1p = &g1[0]; const T* g2p = &g2[0]; \
  for (uint i = 0; i != g1.size(); ++i) *g1p++ OP1 *g2p++; \
  return g1; \
} \
template<class T> grid<T> \
operator OP2(const grid<T>& g1, const grid<T>& g2) \
  { grid<T> tmp = g1; return operator OP1(tmp, g2); }

GGOP(+=, +) GGOP(-=, -) GGOP(*=, *) GGOP(/=, /)

#undef GGOP

// Add/Subtract/Multiply/Divide each element of a vector or grid 'g' to/by
// a fixed scalar value 't'

#define VGSOP1(VG, OP1, OP2) \
template<class T, class S> VG <T>& \
operator OP1(VG <T>& g, const S& t) { \
  T *gp = &g[0]; \
  for (uint i = 0; i != g.size(); ++i) *gp++ OP1 t; \
  return g; \
} \
template<class T, class S> VG <T> \
operator OP2(const VG <T>& g, const S& t) \
  { VG <T> tmp = g; return operator OP1(tmp, t); }

#define VGSOP(OP1, OP2) VGSOP1(vector, OP1, OP2) VGSOP1(grid, OP1, OP2)

VGSOP(+=, +) VGSOP(-=, -) VGSOP(*=, *) VGSOP(/=, /)

#undef VGSOP1
#undef VGSOP

// ==========================================================================
// Mathematical Functions.  Why doesn't 'transform' work for these functions?
// More efficient: function(a, b) --- Less efficient: b = function(a)

// log(grid, grid), log(grid, vector), grid = log(grid), vector = log(vector)
#define VGFN1(FN) \
template <class T> grid <T>& FN(const grid <T>& m1, grid <T>& m2) { \
  m2.resize(m1.rows(), m1.cols()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
} \
template <class T> vector <T>& FN(const grid <T>& m1, vector <T>& m2) { \
  m2.resize(m1.size()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
}

// log(vector, grid), log(vector, vector)
#define VGFN2(VG, FN) \
template <class T> VG <T>& FN(const vector <T>& m1, VG <T>& m2) { \
  m2.resize(m1.size()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
} \
template<class T> VG <T> FN(const VG <T>& m) \
  { VG <T> tmp; return FN(m, tmp); }

#define VGFN(FN) VGFN1(FN) VGFN2(grid, FN) VGFN2(vector, FN)

VGFN(abs) VGFN(exp) VGFN(log) VGFN(log10) VGFN(sin) VGFN(cos) VGFN(tan)
VGFN(asin) VGFN(acos) VGFN(atan) VGFN(sinh) VGFN(cosh) VGFN(tanh)

#undef VGFN1
#undef VGFN2
#undef VGFN

template <class T>
grid<T>& pow(const grid<T>& m1, T p, grid<T>& m2) {
  m2.resize(m1.rows(), m1.cols());
  typename vector<T>::const_iterator m1i = m1.begin();
  typename vector<T>::iterator m2i = m2.begin();
  while (m1i != m1.end()) *m2i++ = std::pow(*m1i++, p);
  return m2;
}

template<class T> grid<T> pow(const grid<T>& m, T p)
  { grid<T> tmp; return pow(m, p, tmp); }

// ==========================================================================
// Vector/Grid transpose.

template<class T> grid<T>& trans(const grid<T>& g1, grid<T>& g2)
{
  if (&g1 == &g2)
    { grid<T> tmp; g2.swap(trans(g1, tmp)); }
  else {
    g2.resize(g1.cols(), g1.rows());
    for (uint i = 0; i != g1.cols(); ++i)
      for (uint j = 0; j != g1.rows(); ++j)
        g2(i, j) = g1(j, i);
  }
  return g2;
}

// Syntax 2: Nicer syntax atrans = trans(a);
template<class T> grid<T> trans(const grid<T>& g)
  { grid<T> tmp; return trans(g, tmp); }

template<class T> grid<T>& trans(const vector<T>& g1, grid<T>& g2)
{
  if (&g1 == &(const vector<T>&)g2)
    g2.resize(1, (uint)g1.size());
  else {
    g2.resize(1, (uint)g1.size());
    std::copy(g1.begin(), g1.end(), g2.begin());
  }
  return g2;
}

// Syntax 2: Nicer syntax atrans = trans(a);
template<class T> grid<T> trans(const vector<T>& g)
  { grid<T> tmp; return trans(g, tmp); }

// ==========================================================================
// Perform a linear transformation on the values of a grid so that the
// values range from 'val1' to 'val2'

template<class T>
void grid<T>::normalize(T val1, T val2)
{
  if (size() == 0) return;

  T gridmin = min(), gridmax = max();
  if (gridmax > gridmin) {
    T* thisp = &front();
    T vscale = val2 - val1, gscale = gridmax - gridmin;
    T translate = gridmax * val1 - gridmin * val2;
    for (uint i = 0; i != size(); ++i) {
      *thisp = (*thisp * vscale + translate) / gscale;
      thisp++;
    }
  }
  else
    fill(val1 / 2 + val2 / 2);

} // grid::normalize

// ==========================================================================
// quicksort the rows of a grid according to the values in
// a specified column.

template<class T>
void sortrows(grid<T>& a, uint col=0, uint left=INT_MIN, uint right=INT_MIN)
//void sortrows(grid<T>& a, uint col=0, uint left=-1u, uint right=-1u)
{
  if (left == -1u && right == -1u)
    { left = 0; right = a.rows() - 1; }

  uint i = left, j = right;
  T midval = a((left + right) / 2, col);
  do {
    while (a(i, col) < midval && i < right) ++i;
    while (midval < a(j, col) && j > left) --j;
    if (i <= j) {
      for (uint k = 0; k < a.cols(); ++k)
        swap(a(i, k), a(j, k));
      ++i; --j;
    }
  } while (i <= j);
  if (left < j) sortrows(a, col, left, j);
  if (i < right) sortrows(a, col, i, right);

} // sortrows

// ==========================================================================
// LU Decomposition with partial pivoting and associated functions.
// For use when LAPACK is not available.  The pivots indexing is
// consistent with LAPACK.

// Can be safely called as LU(a, a), which is efficient but replaces
// 'a' with its LU decomposition.

template<class T>
int LU(const grid<T>& a, grid<T>& g, ivector& pivots,
       int gridtype=gengrid, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("LU: Grid is not square");

  uint n = a.rows();
  pivots.resize(n);
  if (gridtype == posgrid)
    return chol(a, g, ul);

  g = a;
  if (gridtype == trigrid)
    return 0;
  
  int info = 0;
  for (uint i = 0; i < n - 1; ++i)
  {
    uint maxrow = i;
    double maxval = std::abs(g(i, i));
    for (uint j = i + 1; j < n; ++j)
      if (std::abs(g(j, i)) > maxval) {
        maxval = std::abs(g(j, i));
        maxrow = j;
      }
    if (maxval == 0) {
      info = 1;
      pivots[i] = 1;
    } else {
      pivots[i] = maxrow + 1;
      if (i != maxrow)
        for (uint j = 0; j < n; ++j)
          std::swap(g(i, j), g(maxrow, j));
      for (uint j = i + 1; j < n; ++j) {
        g(j, i) /= g(i, i);
        for (uint k = i + 1; k < n; ++k)
          g(j, k) -= g(j, i) * g(i, k);
      }
    }
  }
  if (g(n - 1, n - 1) == T(0)) {
    info = 1;
    pivots[n - 1] = 1;
  } else
    pivots[n - 1] = n;
  return info;

} // LU

template<class T>
int LUsolve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
    const ivector& pivots, int gridtype=gengrid, char ul='U',
            char tr='N', char dg='N')
{
  if (gridtype == trigrid)
    return trsolve(a, b, sol, ul, tr, dg);

  if (gridtype == posgrid)
    return cholsolve(a, b, sol, ul);

  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("LUsolve: Grid is not square");
  if (a.rows() != pivots.size())
    throw typename grid<T>::grid_error("LUsolve: Pivots vector is the wrong size");
  if (a.rows() != b.rows())
    throw typename grid<T>::grid_error("LUsolve: Grids are not conformable");

  if (&sol == &a || &sol == &b) {
    grid<T> tmp;
    LUsolve(a, b, tmp, pivots, gridtype, ul, tr, dg);
    sol << tmp;
    return 0;
  }

  if (&a == &b) {
    sol = grid<T>("I", a.rows());
    return 0;
  }

  sol = b;
  for (uint i = 0; i < pivots.size(); ++i)
    if ((uint)pivots[i] != i + 1)
     for (uint j = 0; j != sol.cols(); ++j)
        std::swap(sol(i, j), sol(pivots[i] - 1, j));

  grid<T> y(sol.rows(), sol.cols());
  for (uint j = 0; j != sol.cols(); ++j) {
    for (uint i = 0; i != sol.rows(); ++i) {
      y(i, j) = sol(i, j);
      for (uint k = 0; k < i; ++k)
        y(i, j) -= a(i, k) * y(k, j);
    }
  }

  for (uint j = 0; j != sol.cols(); ++j) {
    for (uint i = sol.rows(); i > 0; --i) {
      sol(i - 1, j) = y(i - 1, j) / a(i - 1, i - 1);
      for (uint k = i; k != sol.rows(); ++k)
        sol(i - 1, j) -= a(i - 1, k) * sol(k, j) / a(i - 1, i - 1);
    }
  }
  return 0;
}

template<class T>
int solve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
          int gridtype=gengrid, char ul='U', char tr='N', char dg='N')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("solve: Grid is not square");
  if (a.rows() != b.rows())
    throw typename grid<T>::grid_error("solve: Grids are not conformable");

  ivector pivots; int info = 0;
  if (gridtype == gengrid || gridtype == symgrid || gridtype == posgrid) {
    grid<T> LUthis;
    info = LU(a, LUthis, pivots, gridtype, ul);
    if (info == 0)
      info = LUsolve(LUthis, b, sol, pivots, gridtype, ul, tr, dg);
  } else if (gridtype == trigrid) {
    info = LUsolve(a, b, sol, pivots, gridtype, ul, tr, dg);
  }
  return info;
}

template<class T>
T det(const grid<T>& a, int gridtype=gengrid)
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("det: Grid is not square");
  if (gridtype == trigrid)
    return detLU(a, gridtype);
  if (gridtype == gengrid || gridtype == symgrid || gridtype == posgrid) {
    grid<T> lu; ivector pivots;
    uint info = LU(a, lu, pivots, gridtype);
    if (info != 0) {
      std::cerr << "Warning: LU routine returned info=" << info
           << " with gridtype=" << gridtype << " in det()\n";
      return 0;
    }
    return detLU(lu, gridtype, pivots);
  } else
    throw typename grid<T>::grid_error("det: Unknown grid type");
  return 0;
}

template<class T>
T detLU(const grid<T>& lu, int gridtype=gengrid, const ivector pivots=ivector())
{
  if (lu.rows() != lu.cols())
    throw typename grid<T>::grid_error("detLU: Grid is not square");
  T res = 1;
  if (gridtype == gengrid || gridtype == symgrid) {
#ifdef USE_LAPACK
    if (gridtype == gengrid) {
#else
    if (gridtype == gengrid || gridtype == symgrid) {
#endif
      int mult = 1;
      for (uint i = 0; i < lu.rows(); ++i) {
        res *= lu(i, i);
        if ((uint)pivots[i] != 1+i) mult = -mult;
      }
      res *= mult;
    } else // THIS DOES NOT WORK!!!
      for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
  } else if (gridtype == posgrid) {
    for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
    res *= res;
  } else if (gridtype == trigrid) {
    for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
  } else
    throw typename grid<T>::grid_error("detLU: Unknown grid type");
  return res;
}

template<class T>
int inv(const grid<T>& a, grid<T>& x, int gridtype=gengrid,
        char ul='U', bool fill=true, char dg='N')
{
  grid<T> lu;
  ivector pivots;
  if (LU(a, lu, pivots, gridtype, ul))
    throw typename grid<T>::grid_error("inv: Grid is singular");
  return LUsolve(lu, grid<T>("I", a.rows()), x, pivots, gridtype, ul);
}

template<class T>
grid<T> inv(const grid<T>& a, int gridtype=gengrid, char ul='U',
                bool fill=true, char dg='N')
  { grid<T> tmp; inv(a, tmp, gridtype, ul, fill, dg); return tmp; }//Added fill param to fix conversion warning from passing int gridtype to char ul. -SJZ

template<class T> grid<T>&
backsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='U')
{
  x = b;
  if (ul == 'U') {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = x.rows() - 1 - i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = 0; k != curr; ++k)
          x(k, j) -= a(k, curr) * x(curr, j);
      }
    }
  }
  else {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = x.rows() - 1 - i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = 0; k != curr; ++k)
          x(k, j) -= a(curr, k) * x(curr, j);
      }
    }
  }
  return x;
}

template<class T> grid<T>&
forwardsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='L')
{
  x = b;
  if (ul == 'L') {
    for (uint i = 0; i < x.rows(); ++i) {
      if (a(i, i) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(i, j) /= a(i, i);
        for (uint k = i + 1; k < x.rows(); ++k)
          x(k, j) -= a(k, i) * x(i, j);
      }
    }
  }
  else {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = i + 1; k < x.rows(); ++k)
          x(k, j) -= a(i, k) * x(i, j);
      }
    }
  }
  return x;
}

// ==========================================================================
// Cholesky Decomposition Functions for use when LAPACK is not available.
// Can be called as chol(a, a), which is efficient but replaces 'a' with
// its Cholesky decomposition.  'ul' can be either 'U' or 'L' to
// indicate whether the upper or lower half is computed.

template<class T>
int chol(const grid<T>& a, grid<T>& chol, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("chol: Grid is not square");
  uint n = a.rows();
  chol = a;
  //if (ul == 'U') {		//ALP for UNIX
  //  for (uint i = 0; i != n - 1; ++i)
  //    memset(&chol[0] + i * n + i + 1, 0, sizeof(T) * (n - i - 1));
  //} else {
  //  for (uint i = 1; i != n; ++i)
  //    memset(&chol[0] + i * n, 0, sizeof(T) * i);
  //}
  if (ul == 'U') {
    for (uint i = 0; i != n; ++i)
    {
      for (uint j = 0; j != i; ++j)
        chol(i, i) -= chol(j, i) * chol(j, i);
      if (chol(i, i) == T(0))
        return i;
      chol(i, i) = std::sqrt((T)chol(i, i));
      for (uint j = i + 1; j != n; ++j) {
        for (uint k = 0; k != i; ++k)
          chol(i, j) -= chol(k, i) * chol(k, j);
        chol(i, j) /= chol(i, i);
      }
    }
  } else {
    for (uint i = 0; i != n; ++i)
    {
      for (uint j = 0; j != i; ++j)
        chol(i, i) -= chol(i, j) * chol(i, j);
      if (chol(i, i) == T(0))
        return i;
      chol(i, i) = std::sqrt((T)chol(i, i));
      for (uint j = i + 1; j != n; ++j) {
        for (uint k = 0; k != i; ++k)
          chol(j, i) -= chol(i, k) * chol(j, k);
        chol(j, i) /= chol(i, i);
      }
    }
  }
  return 0;

} // chol

template<class T>
int cholsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("cholsolve: Grid is not square");
  grid<T> y = b;
  forwardsolve(a, b, y, ul);
  backsolve(a, y, x, ul);
  return 0;

} // cholsolve

template<class T>
int cholupdate(grid<T>& a, vector<T>& x)
  { throw typename grid<T>::grid_error("cholupdate: Calling non-Lapack stub"); return 0; }

template<class T>
int trsolve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
            char ul, char tr, char dg)
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("trsolve: Grid is not square");
  if (ul == 'U')
    forwardsolve(a, b, sol, ul);
  else
    backsolve(a, b, sol, ul);
  return 0;

} // trsolve

//============================================================================
// Various vector and matrix norms

template<class T>
double l1vecnorm(const grid<T>& a) {
  double sum = 0; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    sum += std::abs(*tp++);
  return sum;
}

template<class T>
double l2vecnorm(const grid<T>& a) {
  double sum = 0, absterm; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    { absterm = std::abs(*tp++); sum += absterm * absterm; }
  return std::sqrt(sum);
}

template<class T>
double linfvecnorm(const grid<T>& a) {
  double max = 0, absterm; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    { if ((absterm = std::abs(*tp++)) > max) max = absterm; }
  return max;
}

template<class T>
double l1matnorm(const grid<T>& a) {
  double max = 0;
  for (uint j = 0; j < a.cols(); ++j) {
    double sum = 0; const T* tp = &a(0, j);
    for (uint i = 0; i < a.rows(); ++i)
      sum += std::abs(*tp++);
    if (sum > max) max = sum;
  }
  return max;
}

template<class T>
double l2matnorm(const grid<T>& a, int gridtype) {
  grid<T> evecs; dgrid evals;
  eigenvectors(a, evecs, evals, gridtype);
  return linfvecnorm(evals);
}

template<class T>
double linfmatnorm(const grid<T>& a) {
  double max = 0;
  for (uint i = 0; i < a.rows(); ++i) {
    double sum = 0;
    for (uint j = 0; j < a.cols(); ++j)
      sum += std::abs(a(i, j));
    if (sum > max) max = sum;
  }
  return max;
}

template<class T>
double frobnorm(const grid<T>& a) { return l2vecnorm(a); }

// ==========================================================================
// Matrix Multiplication Functions.
//
// Multiply two grids and/or vectors together using matrix multiplication
// as is done with '%*%' in R.  Ideally, code that uses this class should
// be compiled using -DUSE_LAPACK and linked with BLAS/LAPACK, but these
// implementations will allow the code to work regardless.

template<class T> grid<T>&
matmult(const grid<T>& a, const grid<T>& b, grid<T>& c,
        char ta='N', char tb='N')
{
  if (a.cols() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == &a || &c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows(), b.cols()); c.fill(0);
    for (uint i = 0; i < a.rows(); ++i)
      for (uint j = 0; j < b.cols(); ++j) {
        const T *ap = &a(i, 0), *bp = &b(0, j);
        T *cp = &c(i, j);
        for (uint k = 0; k < a.cols(); ++k)
          { *cp += *ap * *bp++; ap += a.rows(); }
      }
  }
  return c;
}

template<class T> grid<T>&
matmult(const grid<T>& a, const vector<T>& b, grid<T>& c, char ta='N', char tb='N')
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == &a || (const vector<T>*)&c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows(), 1); c.fill(0);
    const T *ap = &a[0], *bp = &b[0];
    for (uint i = 0; i < a.cols(); ++i) {
      T* cp = &c[0];
      for (uint j = 0; j < a.rows(); ++j) *cp++ += *ap++ * *bp;
      ++bp;
    }
  }
  return c;
}

template<class T> grid<T>&
matmult(const vector<T>& a, const grid<T>& b, grid<T>& c, char ta='N', char tb='N')
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if ((const vector<T>*)&c == &a || &c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(1, b.cols()); c.fill(0);
    const T *bp = &b[0]; T* cp = &c[0];
    for (uint j = 0; j < b.cols(); ++j) {
      const T *ap = &a[0];
      for (uint i = 0; i < b.rows(); ++i) *cp += *ap++ * *bp++;
      ++cp;
    }
  }
  return c;
}

template<class T> vector<T>&
matmult(const grid<T>& a, const vector<T>& b, vector<T>& c, char ta='N', char tb='N')
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == (const vector<T>*)&a || &c == &b)
    { vector<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows()); std::fill(c.begin(), c.end(), 0);
    const T *ap = &a[0], *bp = &b[0];
    for (uint i = 0; i < a.cols(); ++i) {
      T* cp = &c[0];
      for (uint j = 0; j < a.rows(); ++j) *cp++ += *ap++ * *bp;
      ++bp;
    }
  }
  return c;
}

template<class T> vector<T>&
matmult(const vector<T>& a, const grid<T>& b, vector<T>& c, char ta='N', char tb='N')
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (&c == &a || &c == (const vector<T>*)&b)
    { vector<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(b.cols()); std::fill(c.begin(), c.end(), 0);
    T* cp = &c[0];
    const T *bp = &b[0];
    for (uint j = 0; j < b.cols(); ++j) {
      const T* ap = &a[0];
      for (uint i = 0; i < b.rows(); ++i)
        *cp += *ap++ * *bp++;
      ++cp;
    }
  }
  return c;
}

template<class T, class S> grid<T>&
diagmult(const grid<T>& a, const vector<S>& b, grid<T>& c)
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("diagmult: Grids are not conformable");
  if (&c == &a || (const vector<S>*)&c == &b)
    { grid<T> tmp; c.swap(diagmult(a, b, tmp)); }
  else {
    c.resize(a.rows(), a.cols());
    const T* ap = &a[0]; const S* bp = &b[0];
    T* cp = &c[0];
    for (uint i = 0; i < a.size(); ) {
      *cp++ = *ap++ * *bp;
      if (++i % a.rows() == 0) ++bp;
    }
  }
  return c;
}

template<class T, class S> grid<T>&
diagmult(const vector<S>& a, const grid<T>& b, grid<T>& c)
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("diagmult: Grids are not conformable");
  if ((const vector<S>*)&c == &a || &c == &b)
    { grid<T> tmp; c.swap(diagmult(a, b, tmp)); }
  else {
    c.resize(b.rows(), b.cols());
    const S* ap = &a[0]; const T* bp = &b[0];
    T* cp = &c[0];
    for (uint i = 0; i < b.size(); ) {
      *cp++ = *ap++ * *bp++;
      if (++i % b.rows() == 0) ap = &a[0];
    }
  }
  return c;
}

template<class T> grid<T>
matmult(const grid<T>& a, const grid<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T> grid<T>
matmult(const grid<T>& a, const vector<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T> grid<T>
matmult(const vector<T>& a, const grid<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T, class S> grid<T>
diagmult(const grid<T>& a, const vector<S>& b)
  { grid<T> tmp; return diagmult(a, b, tmp); }

template<class T, class S> grid<T>
diagmult(const vector<S>& a, const grid<T>& b)
  { grid<T> tmp; return diagmult(a, b, tmp); }

// ==========================================================================
// Statistical Functions.  Most of these functions accept an argument 'rc'
// that is either "row", "col" or "".  If "row", the next argument(s)
// is the row(s) of the grid that the particular statistical function
// should be applied to.  Similarly, for "col".  If "rc" is empty,
// the function will be applied to all the values in the grid.  If
// an additional argument 'l' is supplied, it is a list (subset)
// of the elements within the row or column that are to be used in
// the computation.

template <class T>
T sum(const grid<T>& g, const string& rc="", uint ind=0)
{
  T total = 0.0;
  const T* curr;
  if (rc == "row") {
    curr = &g(ind, 0);
    for (uint i = 0; i != g.cols(); ++i) {
      total += *curr; curr += g.rows();
    }
  } else if (rc == "col") {
    curr = &g(0, ind);
    total = std::accumulate(curr, curr + g.rows(), T(0));
  } else {
    total = std::accumulate(g.begin(), g.end(), T(0));
  }
  return total;
}

template <class T>
T sum(const grid<T>& g, const string& rc, uint ind, const uivector& l)
{
  T total = 0.0;
  if (rc == "row")
    for (uint i = 0; i < l.size(); ++i)
      total += g(ind, l[i]);
  else if (rc == "col")
    for (uint i = 0; i < l.size(); ++i)
      total += g(l[i], ind);
  else
    for (uint i = 0; i < l.size(); ++i)
      total += g(l[i]);
  return total;

}

template <class T>
T mean(const grid<T>& g, const string& rc="",
 uint ind=0)
{
  T total = 0.0;
  if (rc == "row" && g.cols() != 0)
    total = sum(g, rc, ind) / T(g.cols());
  else if (rc == "col" && g.rows() != 0)
    total = sum(g, rc, ind) / T(g.rows());
  else if (g.size() != 0)
    total = sum(g) / T(g.size());
  return total;
}

template <class T>
T mean(const grid<T>& g, const string& rc, uint ind, const uivector& l)
{
  T total = 0.0;
  if (l.size() != 0) {
    total = sum(g, rc, ind, l);
    total /= l.size();
  }
  return total;
}

template <class T>
grid<T>& mean(const grid<T>& g, grid<T>& m, const string& rc="col")
{
  if (rc == "row") m.resize(g.rows());
  else if (rc == "col") m.resize(g.cols());
  else throw typename grid<T>::grid_error("mean: Need to specify rc (row or col)");

  for (uint i = 0; i != m.size(); ++i) m[i] = mean(g, rc, i);
  return m;
}

template <class T>
T var(const grid<T>& g, const string& rc="", uint ind=0)
{
  T total = 0, m = mean(g, rc, ind);
  const T* curr;
  if (rc == "row" && g.cols() > 1) { 
    curr = &g(ind, 0);
    for (uint i = 0; i != g.cols(); ++i)
      { total += (*curr - m) * (*curr - m); curr += g.rows(); }
    total /= (g.cols() - 1.0);
  } else if (rc == "col" && g.rows() > 1) { 
    curr = &g(0, ind);
    for (uint i = 0; i != g.rows(); ++i)
      { total += (*curr - m) * (*curr - m); curr++; }
    total /= (g.rows() - 1.0);
  } else if (g.size() > 1) {
    curr = &g.front();
    for (uint i = 0; i != g.size(); ++i)
      { total += (*curr - m) * (*curr - m); curr++; }
    total /= (g.size() - 1.0);
  }
  return total;
}

template <class T>
T var(const grid<T>& g, const string& rc, uint ind, const uivector& l)
{
  if (l.size() < 2) return 0.0;
  T total = 0, m = mean(g, rc, ind);
  if (rc == "row")
    for (uint i = 0; i < l.size(); ++i)
      total += (g(ind, l[i]) - m) * (g(ind, l[i]) - m);
  else if (rc == "col")
    for (uint i = 0; i < l.size(); ++i)
      total += (g(l[i], ind) - m) * (g(l[i], ind) - m);
  else
    for (uint i = 0; i < l.size(); ++i)
      total += (g[l[i]] - m) * (g[l[i]] - m);
  total /= (l.size() - 1.0);
  return total;
}

template <class T>
grid<T>& var(const grid<T>& g, grid<T>& v, const string& rc="col")
{
  if (rc == "row") v.resize(g.rows());
  else if (rc == "col") v.resize(g.cols());
  else throw typename grid<T>::grid_error("var: Need to specify rc (row or col)");

  for (uint i = 0; i != v.size(); ++i) v[i] = var(g, rc, i);
  return v;
}

template <class T>
T sd(const grid<T>& g, const string& rc="", uint ind=0)
  { return std::sqrt(var(g, rc, ind)); }

template <class T>
T sd(const grid<T>& g, const string& rc, uint ind, const uivector& l)
  { return std::sqrt(var(g, rc, ind, l)); }

template <class T>
grid<T>& sd(const grid<T>& g, grid<T>& s, const string& rc="col")
{
  if (rc == "row") s.resize(g.rows());
  else if (rc == "col") s.resize(g.cols());
  else throw typename grid<T>::grid_error("sd: Need to specify rc (row or col)");

  for (uint i = 0; i != s.size(); ++i) s[i] = sd(g, rc, i);
  return s;
}

template <class T>
T cov(const grid<T>& g, const string& rc, uint i1, uint i2)
{
  T total = 0.0, m1 = mean(g, rc, i1), m2 = mean(g, rc, i2);
  const T *curr1, *curr2;
  if (rc == "row") {
    curr1 = &g(i1, 0); curr2 = &g(i2, 0);
    for (uint i = 0; i < g.cols(); ++i) {
      total += (*curr1 - m1) * (*curr2 - m2);
      curr1 += g.rows(); curr2 += g.rows();
    }
    total /= (g.cols() - 1.0);
  } else if (rc == "col") {
    curr1 = &g(0, i1); curr2 = &g(0, i2);
    for (uint i = 0; i < g.rows(); ++i)
      total += (*curr1++ - m1) * (*curr2++ - m2);
    total /= (g.rows() - 1.0);
  } else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");
  return total;
}

template <class T>
T cov(const grid<T>& g, const string& rc, uint i1, uint i2, const uivector& l)
{
  T total = 0.0, m1 = mean(g, rc, i1), m2 = mean(g, rc, i2);
  if (rc == "row") {
    for (uint i = 0; i < l.size(); ++i)
      total += (g(i1, l[i]) - m1) * (g(i2, l[i]) - m2);
    total /= (l.size() - 1.0);
  } else if (rc == "col") {
    for (uint i = 0; i < l.size(); ++i)
      total += (g(l[i], i1) - m1) * (g(l[i], i2) - m2);
    total /= (l.size() - 1.0);
  } else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");
  return total;
}

template <class T>
grid<T>& cov(const grid<T>& g, grid<T>& d, const string& rc="col")
{
  int nvars = 0, nobs = 0;
  if (rc == "row")
    { nvars = g.rows(); nobs = g.cols(); matmult(g, trans(g), d); }
  else if (rc == "col")
    { nvars = g.cols(); nobs = g.rows(); matmult(trans(g), g, d); }
  else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");

  grid<T> m(nvars, 1); mean(g, m, rc);
  ((d /= nobs) -= matmult(m, trans(m))) *= (nobs / (nobs - 1.0));
  return d;
}

template <class T>
T cor(const grid<T>& g, const string& rc, uint i1, uint i2)
  { return cov(g, rc, i1, i2) / (sd(g, rc, i1) * sd(g, rc, i2)); }

template <class T>
T cor(const grid<T>& g, const string& rc,
      uint i1, uint i2, const uivector& l)
  { return cor(g, rc, i1, i2, l) / (sd(g, rc, i1, l) * sd(g, rc, i2, l)); }

template <class T>
void cor(const grid<T>& g, grid<T>& d, const string& rc="col")
{
  uint sz = 0;
  if (rc == "row") sz = g.rows();
  else if (rc == "col") sz = g.cols();
  else throw typename grid<T>::grid_error("cor: Need to specify rc (row or col)");
  d.resize(sz, sz);
  for (uint i = 0; i != sz; ++i)
    for (uint j = 0; j != sz; ++j)
      d(i, j) = cor(g, rc, i, j);
}

// ==========================================================================
// Stubs for functions that require BLAS/LAPACK for their implementation.

template<class T>
int eigenvalues(const grid<T>& a, grid<double>& m,
                int gridtype=symgrid, char ul='U')
  { grid<T> e; return eigenvectors(a, e, m, gridtype, ul); }

template<class T>
int eigenvalues(const grid<T>& a, grid<cplx>& m, int gridtype=gengrid,
                char ul='U')
  { grid<T> e; dgrid d; return eigenvectors(a, e, d, gridtype, ul); }

template<class T> int eigenvectors(const grid<T>& a, grid<T>& E, grid<double>& D,
                                   int gridtype=symgrid, char ul='U') {
  throw typename grid<T>::grid_error("eigenvectors: Calling non-Lapack stub");
  return 0;
}

template<class T>
int singularvalues(const grid<T>& a, grid<double>& evals)
  { grid<T> evecs; return singularvectors(a, evecs, evals); }

template<class T>
int singularvectors(const grid<T>& a, grid<T>& E, grid<double>& D) {
  throw typename grid<T>::grid_error("singularvectors: Calling non-Lapack stub");
  return 0;
}

// ==========================================================================
// Specializations to dgrid/zgrid and additional functions for use
// when BLAS/LAPACK libraries are available

#ifdef USE_LAPACK
	template<> dgrid& matmult(const dgrid& a, const dgrid& b, dgrid& c, char ta, char tb);
	template<> dgrid& matmult(const dgrid& a, const dvector& b, dgrid& c, char ta, char tb);
	template<> dgrid& matmult(const dvector& a, const dgrid& b, dgrid& c, char ta, char tb);
	template<> dvector& matmult(const dgrid&, const dvector&, dvector& c, char ta, char tb);
	template<> dvector& matmult(const dvector&, const dgrid&, dvector& c, char ta, char tb);
	template<> zgrid& matmult(const zgrid& a, const zgrid& b, zgrid& c, char ta, char tb);
	template<> zgrid& matmult(const zgrid& a, const zvector& b, zgrid& c, char ta, char tb);
	template<> zgrid& matmult(const zvector& a, const zgrid& b, zgrid& c, char ta, char tb);
	template<> zvector& matmult(const zgrid&, const zvector&, zvector& c, char ta, char tb);
	template<> zvector& matmult(const zvector&, const zgrid&, zvector& c, char ta, char tb);

	template<> int LU(const dgrid&, dgrid&, ivector&, int, char);
	template<> int LUsolve(const dgrid&, const dgrid&, dgrid&, const ivector&,
						   int, char, char, char);
	template<> int inv(const dgrid&, dgrid&, int, char, bool, char);
	template<> dgrid inv(const dgrid&, int, char, bool, char);
	template<> dgrid& backsolve(const dgrid&, const dgrid&, dgrid&, char);
	template<> dgrid& forwardsolve(const dgrid&, const dgrid&, dgrid&, char);
	template<> int chol(const dgrid&, dgrid&, char);
	template<> int cholsolve(const dgrid&, const dgrid&, dgrid&, char);
	template<> int cholupdate(dgrid& a, dvector& x);
	template<> int cholupdate(zgrid& a, zvector& x);
	template<> int trsolve(const dgrid&, const dgrid&, dgrid&, char, char, char);
	template<> int eigenvalues(const dgrid&, dgrid&, int, char ul);
	template<> int eigenvalues(const dgrid&, zgrid&, int, char ul);
	template<> int eigenvectors(const dgrid&, dgrid&, dgrid&,
								int gridtype, char ul);
	template<> int singularvalues(const dgrid&, dgrid&);
	template<> int singularvectors(const dgrid&, dgrid&, dgrid&);

	template<> int LU(const zgrid&, zgrid&, ivector&, int, char);
	template<> int LUsolve(const zgrid&, const zgrid&, zgrid&, const ivector&, int,
	  char, char, char);
	template<> int inv(const zgrid&, zgrid&, int, char, bool, char);
	template<> zgrid inv(const zgrid&, int, char, bool, char);
	template<> int chol(const zgrid&, zgrid&, char);
	template<> int cholsolve(const zgrid&, const zgrid&, zgrid&, char);
	template<> int trsolve(const zgrid&, const zgrid&, zgrid&, char, char, char);

	template<> zgrid& backsolve(const zgrid&, const zgrid&, zgrid&, char);
	template<> zgrid& forwardsolve(const zgrid&, const zgrid&, zgrid&, char);
	template<> int eigenvalues(const zgrid&, dgrid&, int, char ul);
	template<> int eigenvalues(const zgrid&, zgrid&, int, char ul);
	template<> int eigenvectors(const zgrid&, zgrid&, dgrid&, int, char ul);
	template<> int singularvalues(const zgrid&, dgrid&);
	template<> int singularvectors(const zgrid&, zgrid&, dgrid&);

	// ==========================================================================
	// BLAS routines

	dgrid& axpy(const dgrid& x, double alpha, dgrid& y);
	dgrid& gemm(const dgrid& a, const dgrid& b, dgrid& c, double alpha=1,
				double beta=0, char tra='N', char trb='N');
	dgrid& trmm(const dgrid& a, dgrid& b, double alpha=1, const char side='L',
				const char uplo='U', const char tra='N', const char diag='N');
	dgrid& trsm(const dgrid& a, const dgrid& b, dgrid& x, double alpha,
				char sd='L', char ul='U', char tr='N', char dg='N');

	zgrid& axpy(const zgrid& x, const cplx alpha, zgrid& y);
	zgrid& gemm(const zgrid& a, const zgrid& b, zgrid& c, const cplx& alpha=1,
				const cplx& beta=0, char tra='N', char trb='N');
	zgrid& trmm(const zgrid& a, zgrid& b, const cplx& alpha=1, const char side='L',
				const char uplo='U', const char tra='N', const char diag='N');
	zgrid& trsm(const zgrid& a, const zgrid& b, zgrid& x, const cplx& alpha,
				char sd='L', char ul='U', char tr='N', char dg='N');

	// ==========================================================================
	// LAPACK routines

	// Routines to solve AX=B.  A is replaced by some LU type of
	// decomposition, while B is replaced by the solution X.

	int gels(dgrid& a, dgrid& b, char trans='N');
	int gesv(dgrid& a, dgrid& b, ivector& pivots);
	int sysv(dgrid& a, dgrid& b, ivector& pivots, char ul);
	int posv(dgrid& a, dgrid& b, char ul);

	int gels(zgrid& a, zgrid& b, char trans='N');
	int gesv(zgrid& a, zgrid& b, ivector& pivots);
	int sysv(zgrid& a, zgrid& b, ivector& pivots, char ul);
	int hesv(zgrid& a, zgrid& b, ivector& pivots, char ul);
	int posv(zgrid& a, zgrid& b, char ul);

#endif

// ==========================================================================

#ifdef GRIDTEST
template class grid<unsigned char>;
template class grid<int>;
template class grid<double>;
//template class grid<cplx>;
#endif // GRIDTEST

//#include "nr.h"

// ==========================================================================
// grid.h

