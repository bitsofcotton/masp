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
#include <cmath>
#include <limits>
#include <assert.h>
#include <stdint.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;
lieonn_t lieonn;

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
  lieonnStaticInit();
  if(argc <= 1 || argv[1][1]) goto usage;
  if(m == '-') {
    SimpleMatrix<num_t> L;
    std::cin >> L;
    if(2 < argc) {
      for(int i0 = 2; i0 < argc; i0 ++) {
        vector<SimpleMatrix<num_t> > in;
        if(! loadp2or3<num_t>(in, argv[i0])) exit(- 1);
        for(int i = 0; i < in.size(); i ++) {
          SimpleMatrix<num_t> wm(in[i].rows(), in[i].cols() * 9);
          wm.O();
          for(int j = 0; j < in[i].rows(); j ++)
            wm.row(j) = bitsSlide<num_t>(in[i].row(j), 9);
          in[i] = move(wm);
        }
        const auto out4(maspGet<num_t>(in, L));
        vector<SimpleMatrix<num_t> > out;
        out.resize(in.size(), SimpleMatrix<num_t>(1, out4.size()));
        for(int j = 0; j < out.size(); j ++)
          out[j].row(0) = out4;
        if(! savep2or3<num_t>((string(argv[i0]) + string("-4.ppm")).c_str(), out))
          cerr << "failed to save." << endl;
      }
    } else {
      std::string s;
      std::getline(std::cin, s, '\n');
      while(std::getline(std::cin, s, '\n')) {
        SimpleVector<num_t> out(maspGet<num_t>(s2sv<num_t>(s), L));
        for(int i = 0; i < out.size() - 1; i ++)
          std::cout << out[i] << ", ";
        std::cout << out[out.size() - 1] << std::endl;
      }
    }
  } else if(m == 'i') {
    SimpleMatrix<num_t> L;
    std::cin >> L;
    if(2 < argc) {
      assert(3 < argc);
      const int height(std::atoi(argv[2]));
      for(int i0 = 3; i0 < argc; i0 ++) {
        vector<SimpleMatrix<num_t> > in;
        if(! loadp2or3<num_t>(in, argv[i0])) exit(- 1);
        assert(in[0].rows() == 1);
        in = maspApply<num_t>(in, L, height);
        for(int i = 0; i < in.size(); i ++) {
          SimpleMatrix<num_t> wm(in[i].rows(), in[i].cols() / 9);
          wm.O();
          for(int j = 0; j < in[i].rows(); j ++)
            wm.row(j) = bitsG<num_t, false>(in[i].row(j), - 9);
          in[i] = move(wm);
        }
        if(! savep2or3<num_t>((string(argv[i0]) + string("-i.ppm")).c_str(),
          normalize<num_t>(in)) ) cerr << "failed to save." << endl;
      }
    } else {
      std::string s;
      std::getline(std::cin, s, '\n');
      while(std::getline(std::cin, s, '\n')) {
        SimpleVector<num_t> out(maspApply<num_t>(s2sv<num_t>(s), L));
        for(int i = 0; i < out.size() - 1; i ++)
          std::cout << out[i] << ", ";
        std::cout << out[out.size() - 1] << std::endl;
      }
    }
  } else if(m == '+') {
    if(2 < argc) {
      vector<vector<SimpleMatrix<num_t> > > in;
      in.reserve(argc - 1);
      for(int i0 = 2; i0 < argc; i0 ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i0])) continue;
        for(int i = 0; i < work.size(); i ++) {
          SimpleMatrix<num_t> wm(work[i].rows(), work[i].cols() * 9);
          wm.O();
          for(int j = 0; j < work[i].rows(); j ++)
            wm.row(j) = bitsSlide<num_t>(work[i].row(j), 9);
          work[i] = move(wm);
        }
        in.emplace_back(work);
        assert(work.size()    == in[0].size() &&
               work[0].rows() == in[0][0].rows() &&
               work[0].cols() == in[0][0].cols());
      }
      SimpleMatrix<num_t> mres(maspPut<num_t>(in));
      static bool shown(false);
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
    } else {
      std::vector<SimpleVector<num_t> > work;
      std::string s;
      while(std::getline(std::cin, s, '\n'))
        work.emplace_back(s2sv<num_t>(s));
      if(work[work.size() - 1].size() <= 1) work.resize(work.size() - 1);
      const int loop(num_t(int(2)) * log(num_t(work.size())) / log(num_t(int(2))) );
      cout << maspPut<num_t>(work, loop) << endl;
    }
  } else if(m == 'I') {
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      for(int i = 0; i < work.size(); i ++)
        for(int j = 0; j < work[i].rows(); j ++)
          work[i].row(j) = bitsG<num_t, false>(bitsSlide<num_t>(
            work[i].row(j), 9), - 9);
      if(! savep2or3<num_t>((string(argv[i0]) + string("-id.ppm")).c_str(),
        normalize<num_t>(work)) ) cerr << "failed to save." << endl;
    }
  } else goto usage;
  lieonnStaticDestroy();
  return 0;
 usage:
  lieonnStaticDestroy();
  cerr << "Usage:" << endl;
  cerr << argv[0] << " + in0.ppm ... > L.txt" << endl;
  cerr << argv[0] << " - another0.ppm ... < L.txt" << endl;
  cerr << "ddpmopt(32)?(mp)? p another0.ppm-4.ppm ..." << endl;
  cerr << argv[0] << " i <another0.ppm-height> predg.ppm ... < L.txt" << endl;
  cerr << argv[0] << " + < in.txt > L.txt" << endl;
  cerr << "cat L.txt in.txt | " << argv[0] << " - > masp.csv" << endl;
  cerr << "cat L.txt masp.csv | " << argv[0] << " i > out.txt" << endl;
  return - 1;
}

