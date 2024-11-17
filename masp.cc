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
  assert(1 < argc);
  const auto& m(argv[1][0]);
  assert((m == '+' || m == '-' || m == 'i' || m == '4') && ! argv[1][1]);
  if(m == '-' || m == '4') {
    assert(argc == 3);
    vector<SimpleMatrix<num_t> > in;
    if(! loadp2or3<num_t>(in, argv[2])) exit(- 1);
    SimpleVector<num_t> mi(in.size() * in[0].rows() * in[0].cols());
    mi.O();
    for(int i = 0; i < in.size(); i ++)
      for(int j = 0; j < in[i].rows(); j ++)
         mi.setVector(i * in[i].rows() * in[i].cols() + j * in[i].cols(),
           in[i].row(j));
    mi = makeProgramInvariant<num_t>(mi, - num_t(int(1)), true).first;
    SimpleMatrix<num_t> out4(4, mi.size());
    for(int i = 0; i < mi.size() - 1; i ++) {
      SimpleMatrix<num_t> L;
      std::cin >> L;
      out4.setCol(i, L * mi);
    }
    if(m == '4') {
      out4 = out4 * out4.transpose();
      out4 = out4.SVD() * out4 * out4.transpose().SVD().transpose();
      vector<SimpleMatrix<num_t> > out;
      out.resize(1, SimpleMatrix<num_t>(1, out4.rows()));
      for(int i = 0; i < out[0].cols(); i ++)
        out[0](0, i) = out4(i, i);
      if(! savep2or3<num_t>((string(argv[2]) + string("-4.ppm")).c_str(), normalize<num_t>(out)) )
        cerr << "failed to save." << endl;
      return 0;
    }
    for(int i = 0; i < out4.rows(); i ++)
      out4.row(i) = revertProgramInvariant<num_t>(make_pair(
        normalize<num_t>(out4.row(i) ), num_t(int(1))), true);
    vector<SimpleMatrix<num_t> > oimg;
    oimg.resize(in.size());
    for(int i = 0; i < oimg.size(); i ++) {
      oimg[i].resize(in[i].rows() * 2, in[i].cols() * 2);
      oimg[i].O();
      for(int j = 0; j < in[i].rows(); j ++)
        for(int k = 0; k < 4; k ++)
          oimg[i].row(j + (k & 2 ? 0 : in[i].rows())).setVector(
            k & 1 ? 0 : in[i].cols(),
              out4.row(k).subVector(i * in[i].rows() * in[i].cols() +
                j * in[i].cols(), in[i].cols() ) );
    }
    if(! savep2or3<num_t>((string(argv[2]) + string("-i4.ppm")).c_str(), oimg) )
      cerr << "failed to save." << endl;
  } else if(m == 'i') {
    assert(argc == 3);
    vector<SimpleMatrix<num_t> > in;
    if(! loadp2or3<num_t>(in, argv[2])) exit(- 1);
    vector<SimpleMatrix<num_t> > out;
    out.resize(in.size(),
      SimpleMatrix<num_t>(in[0].rows() / 2, in[0].cols() / 2).O());
    for(int i = 0; i < out.size(); i ++)
      for(int j = 0; j < out[0].rows(); j ++)
        for(int k = 0; k <  out[0].cols(); k ++) {
          SimpleMatrix<num_t> L;
          std::cin >> L;
          assert(L.rows() == 4);
          out[i](j, k) = num_t(int(0));
          for(int n = 0; n < L.rows(); n ++)
            out[i](j, k) +=
              L(n, 1 + i * out[0].rows() * out[0].cols() + j * out[0].cols() +
                   k) / L(n, 0) * in[i](j + (n & 2 ? 0 : out[i].rows()),
                                        k + (n & 1 ? 0 : out[i].cols()));
        }
    if(! savep2or3<num_t>((string(argv[2]) + string("-i.ppm")).c_str(), out) )
      cerr << "failed to save." << endl;
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
    const int d(num_t(int(2)) / sqrt(log(num_t(int(2)))) * sqrt(num_t(int(in.size())) * log(num_t(int(in.size()))) ) );
    assert(in[0].size() * in[0][0].rows() * in[0][0].cols() <= in.size());
    assert(in[0].size() < d);
    if(d * 2 <= in[0].size())
      cerr << "we might need crush input but we don't implement, continue..." << endl;
    SimpleMatrix<num_t> L0(in.size(),
      in[0].size() * in[0][0].rows() * in[0][0].cols() + 1);
    L0.O();
    for(int n = 0; n < in.size(); n ++)
      for(int i = 0; i < in[0].size(); i ++)
        for(int j = 0; j < in[0][0].rows(); j ++)
          L0.row(n).setVector(1 +
            i * in[0][0].rows() * in[0][0].cols() +
              j * in[0][0].cols(), in[n][i].row(j));
    static bool shown(false);
    for(int i = 0; i < in[0].size(); i ++)
      for(int j = 0; j < in[0][0].rows(); j ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int k = 0; k < in[0][0].cols(); k ++) {
          vector<SimpleVector<num_t> > res;
          vector<pair<num_t, int> > sute;
          SimpleMatrix<num_t> mres;
          SimpleMatrix<num_t> L;
          SimpleMatrix<num_t> Q;
          res.reserve(4);
          cerr << i << " / " << in[0].size() << ", " << j << " / ";
          cerr << in[0][0].rows() << ", " << k << " / " << in[0][0].cols();
          cerr << endl;
          L.resize(in.size(), in[0].size() * in[0][0].rows() * in[0][0].cols() + 2);
          for(int n = 0; n < L0.rows(); n ++) {
            auto L0n(L0.row(n));
            L0n[0] = in[n][i](j, k);
            L0n[1 + i * in[n][0].rows() * in[n][0].cols() + j * in[n][0].cols() + k] =
              num_t(int(0));
            L.row(n) =
              makeProgramInvariant<num_t>(L0n, - num_t(int(1)), true).first;
          }
          Q = L.QR();
          for(int ii = 0; ii < Q.rows() - 1; ii ++) {
            const auto orth(linearInvariant<num_t>(L));
            const auto n2(orth.dot(orth));
            if(Q.rows() - 5 <= ii) {
              res.emplace_back(orth);
              if(res[res.size() - 1][0] == num_t(int(0)))
                res[res.size() - 1] *= num_t(int(0));
              else
                res[res.size() - 1] /= - res[res.size() - 1][0];
              res[res.size() - 1] = res[res.size() - 1].subVector(1,
                res[res.size() - 1].size() - 1);
            }
            Q.zeroFix(L, sute);
            for(int jj = 0; jj < Q.cols(); jj ++)
              Q.setCol(jj, Q.col(jj) - orth * Q.col(jj).dot(orth) / n2);
          }
          mres.resize(res.size(), res[0].size());
          mres.entity = move(res);
          for(int i = 0; i < mres.rows(); i ++)
            for(int j = 0; j < mres.cols(); j ++)
              if(! isfinite(mres(i, j)) ) {
                if(! shown) {
                  shown = true;
                  cerr << "failed isfinite(mres) replacing with 0." << endl;
                }
                mres(i, j) = num_t(int(0));
              }
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            cout << mres << endl;
          }
        }
      }
  }
  return 0;
}

