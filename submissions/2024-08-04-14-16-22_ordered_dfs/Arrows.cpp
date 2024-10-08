#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
#include <optional>
#ifdef _MSC_VER
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
#include <omp.h>
/* g++ functions */
int __builtin_clz(unsigned int n) { unsigned long index; _BitScanReverse(&index, n); return 31 - index; }
int __builtin_ctz(unsigned int n) { unsigned long index; _BitScanForward(&index, n); return index; }
namespace std { inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); } }
int __builtin_popcount(int bits) {
    bits = (bits & 0x55555555) + (bits >> 1 & 0x55555555);
    bits = (bits & 0x33333333) + (bits >> 2 & 0x33333333);
    bits = (bits & 0x0f0f0f0f) + (bits >> 4 & 0x0f0f0f0f);
    bits = (bits & 0x00ff00ff) + (bits >> 8 & 0x00ff00ff);
    return (bits & 0x0000ffff) + (bits >> 16 & 0x0000ffff);
}
/* enable __uint128_t in MSVC */
//#include <boost/multiprecision/cpp_int.hpp>
//using __uint128_t = boost::multiprecision::uint128_t;
#endif

/** compro io **/
namespace aux {
    template<typename T, unsigned N, unsigned L> struct tp { static void output(std::ostream& os, const T& v) { os << std::get<N>(v) << ", "; tp<T, N + 1, L>::output(os, v); } };
    template<typename T, unsigned N> struct tp<T, N, N> { static void output(std::ostream& os, const T& v) { os << std::get<N>(v); } };
}
template<typename... Ts> std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) { os << '['; aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t); return os << ']'; } // tuple out
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x); // container out (fwd decl)
template<class S, class T> std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) { return os << "[" << p.first << ", " << p.second << "]"; } // pair out
template<class S, class T> std::istream& operator>>(std::istream& is, std::pair<S, T>& p) { return is >> p.first >> p.second; } // pair in
std::ostream& operator<<(std::ostream& os, const std::vector<bool>::reference& v) { os << (v ? '1' : '0'); return os; } // bool (vector) out
std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) { bool f = true; os << "["; for (const auto& x : v) { os << (f ? "" : ", ") << x; f = false; } os << "]"; return os; } // vector<bool> out
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) { bool f = true; os << "["; for (auto& y : x) { os << (f ? "" : ", ") << y; f = false; } return os << "]"; } // container out
template<class T, class = decltype(std::begin(std::declval<T&>())), class = typename std::enable_if<!std::is_same<T, std::string>::value>::type> std::istream& operator>>(std::istream& is, T& a) { for (auto& x : a) is >> x; return is; } // container in
template<typename T> auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) { out << t.stringify(); return out; } // struct (has stringify() func) out
/** io setup **/
struct IOSetup { IOSetup(bool f) { if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); } std::cout << std::fixed << std::setprecision(15); } }
iosetup(true); // set false when solving interective problems
/** string formatter **/
template<typename... Ts> std::string format(const std::string& f, Ts... t) { size_t l = std::snprintf(nullptr, 0, f.c_str(), t...); std::vector<char> b(l + 1); std::snprintf(&b[0], l + 1, f.c_str(), t...); return std::string(&b[0], &b[0] + l); }
/** dump **/
#ifdef _MSC_VER
#define ENABLE_DUMP
#endif
#ifdef ENABLE_DUMP
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
#else
#define dump(...) void(0);
#endif
/** timer **/
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 2.8e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 2.8e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};
/** rand **/
struct Xorshift {
    Xorshift() {}
    Xorshift(uint64_t seed) { reseed(seed); }
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next_u64(); }
    inline uint64_t next_u64() { x ^= x << 7; return x ^= x >> 9; }
    inline uint32_t next_u32() { return next_u64() >> 32; }
    inline uint32_t next_u32(uint32_t mod) { return ((uint64_t)next_u32() * mod) >> 32; }
    inline uint32_t next_u32(uint32_t l, uint32_t r) { return l + next_u32(r - l + 1); }
    inline double next_double() { return next_u32() * e; }
    inline double next_double(double c) { return next_double() * c; }
    inline double next_double(double l, double r) { return next_double(r - l) + l; }
private:
    static constexpr uint32_t M = UINT_MAX;
    static constexpr double e = 1.0 / M;
    uint64_t x = 88172645463325252LL;
};
/** shuffle **/
template<typename T> void shuffle_vector(std::vector<T>& v, Xorshift& rnd) { int n = v.size(); for (int i = n - 1; i >= 1; i--) { auto r = rnd.next_u32(i); std::swap(v[i], v[r]); } }
/** split **/
std::vector<std::string> split(const std::string& str, const std::string& delim) {
    std::vector<std::string> res;
    std::string buf;
    for (const auto& c : str) {
        if (delim.find(c) != std::string::npos) {
            if (!buf.empty()) res.push_back(buf);
            buf.clear();
        }
        else buf += c;
    }
    if (!buf.empty()) res.push_back(buf);
    return res;
}
std::string join(const std::string& delim, const std::vector<std::string>& elems) {
    if (elems.empty()) return "";
    std::string res = elems[0];
    for (int i = 1; i < (int)elems.size(); i++) {
        res += delim + elems[i];
    }
    return res;
}
/** misc **/
template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) { std::fill((T*)array, (T*)(array + N), val); } // fill array
template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

#if 1
inline double get_temp(double stemp, double etemp, double t, double T) {
    return etemp + (stemp - etemp) * (T - t) / T;
};
#else
inline double get_temp(double stemp, double etemp, double t, double T) {
    return stemp * pow(etemp / stemp, t / T);
};
#endif

struct LogTable {
    static constexpr int M = 65536;
    static constexpr int mask = M - 1;
    double l[M];
    LogTable() : l() {
        unsigned long long x = 88172645463325252ULL;
        double log_u64max = log(2) * 64;
        for (int i = 0; i < M; i++) {
            x = x ^ (x << 7);
            x = x ^ (x >> 9);
            l[i] = log(double(x)) - log_u64max;
        }
    }
    inline double operator[](int i) const { return l[i & mask]; }
} log_table;



constexpr int dy[] = { -1, -1, 0, 1, 1, 1, 0, -1 };
constexpr int dx[] = { 0, 1, 1, 1, 0, -1, -1, -1 };

struct Input {

    const int N;
    const std::vector<std::vector<int>> arrows;
    const std::vector<std::vector<int>> mults;

private:
    Input(
        const int N_,
        const std::vector<std::vector<int>>& arrows_,
        const std::vector<std::vector<int>>& mults_
    ) : N(N_), arrows(arrows_), mults(mults_) {}

public:
    static Input load(std::istream& in) {
        int N;
        in >> N;
        assert(8 <= N && N <= 30);
        auto arrows = make_vector(0, N, N);
        auto mults = make_vector(0, N, N);
        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                in >> arrows[y][x] >> mults[y][x];
            }
        }
        return Input(N, arrows, mults);
    }

    static Input load(const int seed) {
        std::ifstream ifs(format("../../tester/in/%d.in", seed));
        return load(ifs);
    }

};

namespace SimpleDFS {

    struct State {

        Timer& timer;

        int N;
        int arrows[32][32];
        int mults[32][32];

        bool used[32][32];
        int ys[900], xs[900];
        int id = 0;
        int score = 0;

        int best_ys[900], best_xs[900];
        int best_id = 0;
        int best_score = 0;

        size_t dfs_call = 0;


        State(const Input& input, Timer& timer_) : timer(timer_) {
            N = input.N;
            for (int y = 0; y < 32; y++) {
                for (int x = 0; x < 32; x++) {
                    arrows[y][x] = mults[y][x] = -1;
                    used[y][x] = false;
                }
            }
            for (int y = 0; y < N; y++) {
                for (int x = 0; x < N; x++) {
                    arrows[y + 1][x + 1] = input.arrows[y][x];
                    mults[y + 1][x + 1] = input.mults[y][x];
                }
            }
        }

        void dfs(int y, int x) {
            if (timer.elapsed_ms() > 9500) return;
            ys[id] = y;
            xs[id] = x;
            id++;
            used[y][x] = true;
            score += id * mults[y][x];
            if (chmax(best_score, score)) {
                dump(id, y, x, best_score);
                std::memcpy(best_ys, ys, sizeof(int) * id);
                std::memcpy(best_xs, xs, sizeof(int) * id);
                best_id = id;
            }
            const auto dir = arrows[y][x];
            for (int i = 1; i < N; i++) {
                int ny = y + dy[dir] * i, nx = x + dx[dir] * i;
                if (arrows[ny][nx] == -1) break;
                if (used[ny][nx]) continue;
                dfs(ny, nx);
            }
            score -= id * mults[y][x];
            used[y][x] = false;
            id--;
        }

    };

    std::vector<std::string> solve(const Input& input) {

        Timer timer;

        std::vector<std::string> moves;

        State state(input, timer);
        state.dfs(input.N / 2 + 1, input.N / 2 + 1);

        for (int i = 0; i < state.best_id; i++) {
            moves.push_back(std::to_string(state.best_ys[i] - 1) + " " + std::to_string(state.best_xs[i] - 1));
        }

        return moves;
    }

}

namespace OrderedDFS {

    // mult=1,2,3,5 �̏��Ō���

    struct State {

        Timer& timer;

        int N;
        int arrows[32][32];
        int mults[32][32];

        bool used[32][32];
        int ys[900], xs[900];
        int id = 0;
        int score = 0;

        int best_ys[900], best_xs[900];
        int best_id = 0;
        int best_score = 0;

        size_t dfs_call = 0;


        State(const Input& input, Timer& timer_) : timer(timer_) {
            N = input.N;
            for (int y = 0; y < 32; y++) {
                for (int x = 0; x < 32; x++) {
                    arrows[y][x] = mults[y][x] = -1;
                    used[y][x] = false;
                }
            }
            for (int y = 0; y < N; y++) {
                for (int x = 0; x < N; x++) {
                    arrows[y + 1][x + 1] = input.arrows[y][x];
                    mults[y + 1][x + 1] = input.mults[y][x];
                }
            }
        }

        void dfs(int y, int x) {
            dfs_call++;
            if (timer.elapsed_ms() > 9500) return;
            ys[id] = y;
            xs[id] = x;
            id++;
            used[y][x] = true;
            score += id * mults[y][x];
            if (chmax(best_score, score)) {
                dump(id, y, x, best_score);
                std::memcpy(best_ys, ys, sizeof(int) * id);
                std::memcpy(best_xs, xs, sizeof(int) * id);
                best_id = id;
            }
            const auto dir = arrows[y][x];

            for (int i = 1;; i++) {
                int ny = y + dy[dir] * i, nx = x + dx[dir] * i;
                if (arrows[ny][nx] == -1) break;
                if (used[ny][nx] || mults[ny][nx] != 1) continue;
                dfs(ny, nx);
            }
            for (int i = 1;; i++) {
                int ny = y + dy[dir] * i, nx = x + dx[dir] * i;
                if (arrows[ny][nx] == -1) break;
                if (used[ny][nx] || mults[ny][nx] != 2) continue;
                dfs(ny, nx);
            }
            for (int i = 1;; i++) {
                int ny = y + dy[dir] * i, nx = x + dx[dir] * i;
                if (arrows[ny][nx] == -1) break;
                if (used[ny][nx] || mults[ny][nx] != 3) continue;
                dfs(ny, nx);
            }
            for (int i = 1;; i++) {
                int ny = y + dy[dir] * i, nx = x + dx[dir] * i;
                if (arrows[ny][nx] == -1) break;
                if (used[ny][nx] || mults[ny][nx] != 5) continue;
                dfs(ny, nx);
            }

            score -= id * mults[y][x];
            used[y][x] = false;
            id--;
        }

    };

    std::vector<std::string> solve(const Input& input) {

        Timer timer;

        std::vector<std::string> moves;

        State state(input, timer);
        state.dfs(input.N / 2 + 1, input.N / 2 + 1);
        //state.dfs(1, 1);

        dump(timer.elapsed_ms(), state.dfs_call);

        for (int i = 0; i < state.best_id; i++) {
            moves.push_back(std::to_string(state.best_ys[i] - 1) + " " + std::to_string(state.best_xs[i] - 1));
        }

        return moves;
    }

}


int main(int argc, char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    const bool LOCAL_MODE = argc > 1 && std::string(argv[1]) == "local";
    const int seed = 7;

    const auto input = [&]() {
        if (LOCAL_MODE) {
            return Input::load(seed);
        }
        return Input::load(std::cin);
    }();

    auto moves = OrderedDFS::solve(input);

    if (LOCAL_MODE) {
        std::ofstream out(format("../../tester/out/%d.out", seed));
        out << moves.size() << '\n';
        for (std::string m : moves) out << m << '\n';
    }
    else {
        std::ostream& out = std::cout;
        out << moves.size() << '\n';
        for (std::string m : moves) out << m << '\n';
    }

    return 0;
}