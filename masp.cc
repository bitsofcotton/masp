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
  const auto& mm(argv[1][1]); 
  assert(m && (mm == 'a' || mm == 'p'));
  if(m == '-') {
    SimpleMatrix<num_t> L;
    if(mm == 'a') {
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
      for(int i = 0; i < mi.size(); i ++) {
        std::cin >> L;
        out4.setCol(i, L * mi);
      }
      for(int i = 0; i < out4.rows(); i ++) {
        out4.row(i) = revertProgramInvariant<num_t>(make_pair(
          makeProgramInvariant<num_t>(normalize<num_t>(out4.row(i)),
            - num_t(int(1)), true).first, num_t(int(1))), true);
      }
      vector<SimpleMatrix<num_t> > oimg;
      oimg.resize(in.size());
      for(int j0 = 0; j0 < 4; j0 ++) {
        for(int i = 0; i < oimg.size(); i ++) {
          oimg[i].resize(in[i].rows(), in[i].cols());
          for(int j = 0; j < oimg[i].rows(); j ++)
            oimg[i].row(j) =
              out4.row(j0).subVector(i * in[i].rows() * in[i].cols() +
                j * in[i].cols(), in[i].cols() );
        }
        if(! savep2or3<num_t>((string(argv[2]) + string("-") + to_string(j0) + string(".ppm")).c_str(), oimg) )
          cerr << "failed to save." << endl;
      }
    } else if(mm == 'p') {
      assert(argc == 5);
      vector<vector<SimpleMatrix<num_t> > > in;
      in.resize(3);
      for(int i = 2; i < argc; i ++) {
        if(! loadp2or3<num_t>(in[i - 2], argv[i])) return - 1;
        assert(in[i - 2].size() == in[0].size() &&
               in[i - 2][0].rows() == in[0][0].rows() &&
               in[i - 2][0].cols() == in[0][0].cols() );
      }
      SimpleVector<num_t> mi(in.size() * in[0].size() * in[0][0].rows() * in[0][0].cols());
      mi.O();
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].size(); j ++)
          for(int k = 0; k < in[i][j].rows(); k ++)
            mi.setVector(i * in[0].size() * in[0][0].rows() *
              in[0][0].cols() + j * in[0][0].rows() * in[0][0].cols() +
              k * in[0][0].cols(), in[i][j].row(k));
      mi = makeProgramInvariant<num_t>(mi, - num_t(int(1)), true).first;
      SimpleMatrix<num_t> out4(4, mi.size());
      for(int i = 0; i < mi.size(); i ++) {
        std::cin >> L;
        out4.setCol(i, L * mi);
      }
      for(int i = 0; i < out4.rows(); i ++)
        out4.row(i) = revertProgramInvariant<num_t>(make_pair(
          makeProgramInvariant<num_t>(normalize<num_t>(out4.row(i)),
            - num_t(int(1)), true).first, num_t(int(1))), true);
      vector<SimpleMatrix<num_t> > oimg;
      oimg.resize(in[0].size());
      for(int j0 = 0; j0 < out4.rows(); j0 ++) {
        for(int i = 0; i < oimg.size(); i ++) {
          oimg[i].resize(in[0][i].rows(), in[0][i].cols());
          for(int j = 0; j < oimg[i].rows(); j ++)
            oimg[i].row(j) =
              out4.row(j0).subVector(i * in[0][i].rows() * in[0][i].cols() +
                j * in[0][i].cols(), in[0][i].cols() );
        }
        if(! savep2or3<num_t>((string("masp-pred-") + to_string(j0) + string(".ppm")).c_str(), oimg) )
          cerr << "failed to save." << endl;
      }
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
    const int d(num_t(int(2)) / sqrt(log(num_t(int(2)))) * sqrt(num_t(int(in.size())) * log(num_t(int(in.size()))) ) );
    assert(in[0].size() * in[0][0].rows() * in[0][0].cols() <= in.size());
    assert(in[0].size() < d);
    if(d * 2 <= in[0].size())
      cerr << "we might need crush input but we don't implement, continue..." << endl;
    SimpleMatrix<num_t> L0(in.size(), mm == 'a' ?
      in[0].size() * in[0][0].rows() * in[0][0].cols() + 1 :
      in[0].size() * in[0][0].rows() * in[0][0].cols() );
    L0.O();
    for(int n = 0; n < in.size(); n ++)
      for(int i = 0; i < in[0].size(); i ++)
        for(int j = 0; j < in[0][0].rows(); j ++)
          L0.row(n).setVector((mm == 'a' ? 1 : 0) +
            i * in[0][0].rows() * in[0][0].cols() +
              j * in[0][0].cols(), in[n][i].row(j));
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
          if(mm == 'a') {
            L.resize(in.size(), in[0].size() * in[0][0].rows() * in[0][0].cols() + 2);
            for(int n = 0; n < L0.rows(); n ++) {
              auto L0n(L0.row(n));
              L0n[0] = in[n][i](j, k);
              L0n[1 + i * in[n][0].rows() * in[n][0].cols() + j * in[n][0].cols() + k] =
                num_t(int(0));
              L.row(n) =
                makeProgramInvariant<num_t>(L0n, - num_t(int(1)), true).first;
            }
          } else if(mm == 'p') {
            L.resize(in.size() - 3,
              2 + 3 * in[0].size() * in[0][0].rows() * in[0][0].cols());
            for(int n = 2; n < L0.rows() - 1; n ++) {
              SimpleVector<num_t> work(1 + L0.cols() * 3);
              for(int ii = 0; ii < 3; ii ++)
                work.setVector(1 + ii * L0.cols(), L0.row(n - 2 + ii));
              work[0] = in[n + 1][i](j, k);
              L.row(n - 2) = makeProgramInvariant<num_t>(work, - num_t(int(1)),
                true).first;
            }
          }
          Q = L.QR();
          for(int ii = 0; ii < (mm == 'a' ? Q.rows() - 1 : Q.rows()); ii ++) {
            const auto orth(linearInvariant<num_t>(L));
            const auto n2(orth.dot(orth));
            if(Q.rows() - (mm == 'a' ? 5 : 4) <= ii) {
              res.emplace_back(orth);
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

