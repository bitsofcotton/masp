#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <random>
#include <assert.h>

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::make_pair;

#include <stdlib.h>

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  const auto& m(argv[1][0]);
  if(argc <= 1 || argv[1][1]) goto usage;
  if(m == '-') {
    SimpleMatrix<num_t> L;
    std::cin >> L;
    assert(L.rows() == 4);
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) exit(- 1);
      SimpleVector<num_t> mi(in.size() * in[0].rows() * in[0].cols());
      mi.O();
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].rows(); j ++)
           mi.setVector(i * in[i].rows() * in[i].cols() + j * in[i].cols(),
             in[i].row(j));
      const auto out4(L * makeProgramInvariant<num_t>(mi).first);
      vector<SimpleMatrix<num_t> > out;
      out.resize(in.size(), SimpleMatrix<num_t>(1, out4.size()));
      for(int j = 0; j < out.size(); j ++)
        out[j].row(0) = out4;
      if(! savep2or3<num_t>((string(argv[i0]) + string("-4.ppm")).c_str(),
        normalize<num_t>(out)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == 'i') {
    assert(3 < argc);
    const int height(std::atoi(argv[2]));
    SimpleMatrix<num_t> L;
    std::cin >> L;
    assert(L.rows() == 4);
    for(int i0 = 3; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) exit(- 1);
      assert(in[0].rows() == 1 && in[0].cols() == 4);
      vector<SimpleMatrix<num_t> > out;
      out.resize(in.size(), SimpleMatrix<num_t>(height, L.cols() / (in.size() * height)).O());
      for(int i = 0; i < out.size(); i ++)
        for(int j = 0; j < out[0].rows(); j ++)
          for(int n = 0; n < L.rows(); n ++)
            out[i].row(j) += L.row(n).subVector(
              i * out[0].rows() * out[0].cols() + j * out[0].cols(),
                out[0].cols()) * in[i](0, n);
      if(! savep2or3<num_t>((string(argv[i0]) + string("-i.ppm")).c_str(),
        normalize<num_t>(out)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      in.emplace_back(work);
      assert(work.size()    == in[0].size() &&
             work[0].rows() == in[0][0].rows() &&
             work[0].cols() == in[0][0].cols());
    }
    cerr << "All images are on memory now." << endl << flush;
    const int d(num_t(int(2)) / sqrt(log(num_t(int(2)))) * sqrt(num_t(int(in.size())) * log(num_t(int(in.size()))) ) );
    assert(in[0].size() * in[0][0].rows() * in[0][0].cols() <= in.size());
    assert(in[0].size() < d);
    if(d * 2 <= in[0].size())
      cerr << "we might need crush input but we don't implement, continue..." << endl;
    SimpleMatrix<num_t> L(in.size(),
      in[0].size() * in[0][0].rows() * in[0][0].cols() + 1);
    L.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int n = 0; n < in.size(); n ++) {
      SimpleVector<num_t> l(in[0].size() * in[0][0].rows() * in[0][0].cols());
      for(int i = 0; i < in[0].size(); i ++)
        for(int j = 0; j < in[0][0].rows(); j ++)
          l.setVector(
            i * in[0][0].rows() * in[0][0].cols() +
              j * in[0][0].cols(), in[n][i].row(j));
      L.row(n) = makeProgramInvariant<num_t>(l).first;
    }
    cerr << "All images on memory are formed now, try to QR." << endl << flush;
    vector<SimpleVector<num_t> > res;
    vector<pair<num_t, int> > sute;
    auto Q(L.QR());
    res.reserve(Q.rows() - 1);
    for(int ii = 0; ii < Q.rows() - 1; ii ++) {
      cerr << ii << " / " << Q.rows() - 1 << endl << flush;
      const auto orth(Q.zeroFix(L, sute));
      const auto n2(orth.dot(orth));
      if(res.size()) {
        auto radd(orth / sqrt(n2));
        if(abs(res[res.size() - 1].dot(radd) - num_t(int(1))) <=
          SimpleMatrix<num_t>().epsilon()) break;
        res.emplace_back(move(radd));
      } else
        res.emplace_back(orth / sqrt(n2));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int jj = 0; jj < Q.cols(); jj ++) {
        Q.setCol(jj, Q.col(jj) - orth * Q.col(jj).dot(orth) / n2);
      }
    }
    static bool shown(false);
    SimpleMatrix<num_t> mres;
    mres.resize(4, res[0].size());
    const auto mrow(min(int(4), int(res.size())) );
    for(int i = 0; i < mrow; i ++)
      mres.row(i) = move(res[i - mrow + res.size()]);
    for(int i = mrow; i < mres.rows(); i ++)
      for(int j = 0; j < mres.cols(); j ++)
        mres(i, j) = num_t(int(0));
    for(int i = 0; i < mres.rows(); i ++)
      for(int j = 0; j < mres.cols(); j ++)
        if(! isfinite(mres(i, j)) ) {
          if(! shown) {
            shown = true;
            cerr << "failed isfinite(mres) replacing with 0." << endl;
          }
          mres(i, j) = num_t(int(0));
        }
    cout << mres << endl;
  } else goto usage;
  return 0;
 usage:
  cerr << "Usage:" << endl;
  cerr << argv[0] << " + in0.ppm ... > L.txt" << endl;
  cerr << argv[0] << " - another0.ppm ... < L.txt" << endl;
  cerr << "ddpmopt(32)?(mp)? p another0.ppm-4.ppm ..." << endl;
  cerr << argv[0] << " i <another0.ppm-height> predg.ppm ... < L.txt" << endl;
  return - 1;
}

