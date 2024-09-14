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
using std::atoi;
using std::string;
using std::to_string;
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;
using std::istringstream;

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
    std::string s;
    vector<SimpleMatrix<num_t> > L;
    while(true) {
      SimpleMatrix<num_t> l;
      std::cin >> l;
      if(l.rows() != 4) break;
      L.emplace_back(move(l));
    }
    if(mm == 'a') for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) continue;
      SimpleVector<num_t> mi(in.size() * in[0].rows() * in[0].cols());
      assert(L.size() == mi.size());
      mi.O();
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].rows(); j ++)
           mi.setVector(i * in[i].rows() * in[i].cols() + j * in[i].cols(),
             in[i].row(j));
      mi = makeProgramInvariant<num_t>(mi, - num_t(int(1)), true).first;
      SimpleMatrix<num_t> out4(4, L[0].cols());
      for(int i = 0; i < mi.size(); i ++) {
        assert(L[i].cols() == mi.size());
        out4.setCol(i, L[i] * mi);
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
        if(! savep2or3<num_t>((string(argv[i0]) + string("-") + to_string(j0) + string(".ppm")).c_str(), oimg) )
          cerr << "failed to save." << endl;
      }
    } else if(mm == 'p') {
      assert(argc - 2 == 4);
      vector<vector<SimpleMatrix<num_t> > > in;
      in.resize(4);
      for(int i = 2; i < argc; i ++) {
        if(! loadp2or3<num_t>(in[i - 2], argv[i])) return - 1;
        assert(in[i - 2].size() == in[0].size() &&
               in[i - 2][0].rows() == in[0][0].rows() &&
               in[i - 2][0].cols() == in[0][0].cols() );
      }
      SimpleVector<num_t> rawin;
      rawin.resize(in.size() * in[0].size() * in[0][0].rows() * in[0][0].cols());
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].size(); j ++)
          for(int k = 0; k < in[i][j].rows(); k ++)
            rawin.setVector(i * in[0].size() * in[0][0].rows() *
              in[0][0].cols() + j * in[0][0].rows() * in[0][0].cols() +
              k * in[0][0].cols(), in[i][j].row(k));
      const auto mi(makeProgramInvariant<num_t>(rawin, - num_t(int(1)), true).first);
      assert(L.size() == mi.size());
      SimpleMatrix<num_t> out4(4, mi.size());
      for(int i = 0; i < mi.size(); i ++) {
        assert(L[i].cols() == mi.size());
        out4.setCol(i, L[i] * mi);
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
    vector<SimpleVector<num_t> > in;
    vector<vector<SimpleMatrix<num_t> > > cache;
    cache.resize(3);
    for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      SimpleVector<num_t> tmp;
      tmp.resize(mm == 'a' ? work.size() * work[0].rows() * work[0].cols() + 1:
        work.size() * work[0].rows() * work[0].cols() * 3 + 1);
      tmp.O();
      if(mm == 'a') {
        for(int i = 0; i < work.size(); i ++)
          for(int j = 0; j < work[i].rows(); i ++)
            tmp.setVector(i * work[i].rows() * work[i].cols() +
              j * work[i].cols(), work[i].row(j));
        for(int i = 0; i < work.size(); i ++)
          for(int j = 0; j < work[i].rows(); j ++)
            for(int k = 0; k < work[i].cols(); k ++) {
              tmp[0] = work[i](j, k);
              tmp[1 + i * work[i].rows() * work[i].cols() +
                j * work[i].cols() + k] = num_t(int(0));
              in.emplace_back(tmp);
              tmp[1 + i * work[i].rows() * work[i].cols() +
                j * work[i].cols() + k] = work[i](j, k);
            }
        continue;
      } else {
        if(cache[0].size() == 0) goto nxtp;
        for(int m = 0; m < cache.size(); m ++) {
          assert(cache[m].size() == cache[0].size() &&
                 cache[m][0].rows() == cache[m][0].rows() &&
                 cache[m][0].cols() == cache[m][0].cols());
          for(int k = 0; k < cache[m].size(); k ++)
            for(int n = 0; n < cache[m][k].rows(); n ++)
              tmp.setVector(1 + m * cache[m].size() * cache[0][0].rows() * cache[0][0].cols() +
                k * cache[m][k].rows() * cache[m][k].cols() +
                  n * cache[0][0].cols(),
                    cache[m][k].row(n) );
        }
        for(int i = 0; i < work.size(); i ++)
          for(int j = 0; j < work[i].rows(); j ++)
            for(int k = 0; k < work[i].cols(); k ++) {
              tmp[0] = work[i](j, k);
              in.emplace_back(tmp);
            }
      }
     nxtp:
      cache[0] = move(cache[1]);
      cache[1] = move(cache[2]);
      cache[2] = move(work);
    }
    const int d(num_t(int(2)) / sqrt(log(num_t(int(2)))) * sqrt(num_t(int(in.size())) * log(num_t(int(in.size()))) ) );
    cerr << in.size() << ", " << in[0].size() << ", " << d << endl;
    const auto isize(in.size() / ((in[0].size() - 1) / (mm == 'a' ? 1 : 3)));
    assert(in[0].size() <= isize);
    //if(isize < d) {
    if(isize < d)
      // assert(0 && "we need to crush input, none implemented.");
      cerr << "we need to crush input, none implemented. fall through." << endl;
    //} else {
    for(int i0 = 0; i0 < in[0].size() / isize; i0 ++) {
      SimpleMatrix<num_t> work;
      work.resize(isize, in[0].size());
      for(int i = 0; i < work.rows(); i ++)
        work.row(i) = in[i0 + i];
      auto Q(work.QR());
      vector<pair<num_t, int> > sute;
      vector<SimpleVector<num_t> > res;
      res.reserve(4);
      for(int i = 0; i < work.rows(); i ++) {
        const auto orth(linearInvariant<num_t>(work));
        const auto n2(orth.dot(orth));
        if(work.rows() - 4 <= i) {
          res.emplace_back(orth);
          res[res.size() - 1] /= - res[res.size() - 1][0];
          res[res.size() - 1] = res[res.size() - 1].subVector(1,
            res[res.size() - 1].size() - 1);
        }
        Q.zeroFix(work, sute);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int j = 0; j < Q.cols(); j ++) {
          Q.setCol(j, Q.col(j) - orth * Q.col(j).dot(orth) / n2);
        }
      }
      SimpleMatrix<num_t> mres;
      mres.resize(res.size(), res[0].size());
      mres.entity = move(res);
      cerr << "OK" << endl;
      cout << mres << endl;
    //}
    }
  }
  return 0;
}

