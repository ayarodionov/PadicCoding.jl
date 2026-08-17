// Benchmark of the optimized P = 2 p-adic coder (C++ port), for comparison
// against `demos/benchmark.jl`'s "fast P=2 N=40" row.
//
//   g++ -O3 -std=c++17 padic2_bench.cpp -o padic2_bench
//   ./padic2_bench
//
// Uses the same skewed 256-symbol source (weight 1/(1 + i/8)) and message
// sizes (10,000 / 100,000 / 1,000,000) as the Julia benchmark, and reports
// best-of-k encode/decode throughput.

#include "padic2.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using Clock = std::chrono::steady_clock;

double elapsed_seconds(Clock::time_point t0, Clock::time_point t1) {
    return std::chrono::duration<double>(t1 - t0).count();
}

template <typename F>
double bench(F f, int k = 5) {
    f();  // warmup
    double best = 1e18;
    for (int i = 0; i < k; ++i) {
        auto t0 = Clock::now();
        f();
        auto t1 = Clock::now();
        double t = elapsed_seconds(t0, t1);
        if (t < best) best = t;
    }
    return best;
}

std::vector<int> make_message(int len, std::mt19937_64& rng) {
    std::vector<double> w(256);
    for (int i = 0; i < 256; ++i) w[i] = 1.0 / static_cast<double>(1 + i / 8);  // integer i/8, matches Julia's i ÷ 8
    std::vector<double> cw(256);
    double sum = 0;
    for (int i = 0; i < 256; ++i) { sum += w[i]; cw[i] = sum; }
    for (double& c : cw) c /= sum;

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::vector<int> msg(len);
    for (int i = 0; i < len; ++i) {
        double u = unif(rng);
        int sym = static_cast<int>(std::lower_bound(cw.begin(), cw.end(), u) - cw.begin());
        msg[i] = sym + 1;  // symbols are 1-based
    }
    return msg;
}

void report(const char* name, int len, double tenc, double tdec, size_t nbits) {
    double bytes = len / 1e6;
    std::printf("%-18s len %8d | encode %8.3f ms %6.2f MB/s | decode %8.3f ms %6.2f MB/s | %5.3f bits/sym\n",
                name, len, 1e3 * tenc, bytes / tenc, 1e3 * tdec, bytes / tdec,
                double(nbits) / len);
}

int main() {
    std::mt19937_64 rng(42);
    const int full_len = 1'000'000;
    std::vector<int> full_msg = make_message(full_len, rng);

    std::vector<uint64_t> freqs(256, 1);
    for (int a : full_msg) freqs[a - 1] += 1;

    std::printf("source: random skewed bytes, %d symbols\n\n", full_len);

    const int N = 40;
    padic2::PadicModel2 m(freqs, N);

    for (int len : {10'000, 100'000, full_len}) {
        std::vector<int> sub(full_msg.begin(), full_msg.begin() + len);

        auto bits = padic2::padic_encode2(m, sub);
        auto dec = padic2::padic_decode2(m, bits);
        if (dec != sub) {
            std::fprintf(stderr, "roundtrip mismatch at len %d\n", len);
            return 1;
        }

        double tenc = bench([&] { padic2::padic_encode2(m, sub); });
        double tdec = bench([&] { padic2::padic_decode2(m, bits); });
        report(("fast   P=2 N=" + std::to_string(N)).c_str(), len, tenc, tdec, bits.nbits);
    }
    return 0;
}
