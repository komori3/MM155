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
//#ifdef _MSC_VER
#define ENABLE_DUMP
//#endif
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
constexpr int NMAX = 32;

struct Input {

    const int N;
    const std::array<std::array<int, NMAX>, NMAX> arrows;
    const std::array<std::array<int, NMAX>, NMAX> mults;

private:
    Input(
        const int N_,
        const std::array<std::array<int, NMAX>, NMAX>& arrows_,
        const std::array<std::array<int, NMAX>, NMAX>& mults_
    ) : N(N_), arrows(arrows_), mults(mults_) {}

public:
    static Input load(std::istream& in) {
        int N;
        in >> N;
        assert(8 <= N && N <= 30);
        std::array<std::array<int, NMAX>, NMAX> arrows, mults;
        for (int y = 0; y < NMAX; y++) {
            for (int x = 0; x < NMAX; x++) {
                arrows[y][x] = mults[y][x] = -1;
            }
        }
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
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



namespace ReversedDFS {

    struct State {

        const int N;
        const std::array<std::array<int, NMAX>, NMAX> arrows;
        const std::array<std::array<int, NMAX>, NMAX> mults;

        Timer& timer;
        double duration_ms;

        // (y, x) に向かって伸びている矢印が存在するセルのリスト
        std::vector<std::vector<std::vector<std::pair<int, int>>>> cands;

        bool used[32][32];
        int ys[900], xs[900];
        int id = 0;
        int score = 0;
        int mult_sum = 0;

        int best_ys[900], best_xs[900];
        int best_id = 0;
        int best_score = 0;

        size_t dfs_call = 0;


        State(const Input& input, Timer& timer_, double duration_ms_) :
            N(input.N), arrows(input.arrows), mults(input.mults), timer(timer_), duration_ms(duration_ms_)
        {
            for (int y = 0; y < 32; y++) {
                for (int x = 0; x < 32; x++) {
                    used[y][x] = false;
                }
            }
            cands = make_vector(std::vector<std::pair<int, int>>(), N + 1, N + 1);
            for (int y = 1; y <= N; y++) {
                for (int x = 1; x <= N; x++) {
                    for (int d = 0; d < 8; d++) {
                        for (int k = 0;; k++) {
                            int ny = y - dy[d] * k, nx = x - dx[d] * k;
                            if (arrows[ny][nx] == -1) break;
                            if (arrows[ny][nx] == d) {
                                cands[y][x].emplace_back(ny, nx);
                            }
                        }
                    }
                    std::sort(cands[y][x].begin(), cands[y][x].end(), [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                        int ma = mults[a.first][a.second], mb = mults[b.first][b.second];
                        int da = abs(y - a.first) + abs(x - a.second), db = abs(y - b.first) + abs(x - b.second);
                        return ma == mb ? da < db : ma > mb; // mult でかい順
                    });
                }
            }
        }

        void dfs(int y, int x) {
            dfs_call++;
            if (timer.elapsed_ms() > duration_ms) return;

            ys[id] = y;
            xs[id] = x;
            id++;
            used[y][x] = true;
            mult_sum += mults[y][x];
            score += mult_sum;

            if (chmax(best_score, score)) {
                dump(id, y, x, best_score);
                std::memcpy(best_ys, ys, sizeof(int) * id);
                std::memcpy(best_xs, xs, sizeof(int) * id);
                best_id = id;
            }

            for (const auto& [ny, nx] : cands[y][x]) {
                if (used[ny][nx]) continue;
                dfs(ny, nx);
            }

            score -= mult_sum;
            mult_sum -= mults[y][x];
            used[y][x] = false;
            id--;
        }

    };

    std::pair<int, std::vector<std::pair<int, int>>> solve(const Input& input, double duration_ms) {

        Timer timer;

        std::vector<std::pair<int, int>> moves;

        State state(input, timer, duration_ms);
        
        const int N = input.N;
        const int cy = input.N / 2 + 1, cx = cy;
        std::vector<std::tuple<int, int, int, int>> cands;
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                int mult = state.mults[y][x];
                int dist = abs(cy - y) + abs(cx - x);
                cands.emplace_back(-mult, dist, y, x);
            }
        }
        std::sort(cands.begin(), cands.end());
        
        {
            auto [mult, dist, sy, sx] = cands[0];
            dump(cands.size(), state.cands[sy][sx].size());
            state.dfs(sy, sx);
        }

        dump(timer.elapsed_ms(), state.dfs_call);

        for (int i = state.best_id - 1; i >= 0; i--) {
            moves.emplace_back(state.best_ys[i], state.best_xs[i]);
        }

        return { state.best_score, moves };
    }

}

namespace PathMutation {

    struct DFSInfo {
        int begin, end, ty, tx;
        std::vector<std::pair<int, int>> old_segment;
        std::vector<std::pair<int, int>> new_segment;
        bool found = false;
        int new_score = 0;
    };

    struct State {

        const Input input;

        std::vector<std::pair<int, int>> move_seq;
        std::vector<int> mult_seq;
        int score = 0;

        std::array<std::array<bool, NMAX>, NMAX> used;

        State(const Input& input_, const std::vector<std::pair<int, int>>& move_seq_)
            : input(input_), move_seq(move_seq_)
        {
            for (auto& v : used) {
                v.fill(false);
            }
            move_seq.insert(move_seq.begin(), std::make_pair(-1, -1));
            mult_seq.push_back(-1);
            for (int i = 1; i < (int)move_seq.size(); i++) {
                const auto& [y, x] = move_seq[i];
                const int mult = input.mults[y][x];
                mult_seq.push_back(mult);
                score += i * mult;
                used[y][x] = true;
            }
        }

        void dfs(DFSInfo& dinfo, int y, int x) {
            assert(!used[y][x]);

            if (y == dinfo.ty && x == dinfo.tx) { // arrived
                if (dinfo.new_segment != dinfo.old_segment) {
                    int nscore = 0;
                    int id = 0;
                    for (int i = 1; i < dinfo.begin; i++) {
                        id++;
                        auto [y2, x2] = move_seq[i];
                        nscore += id * input.mults[y2][x2];
                    }
                    for (const auto& [y2, x2] : dinfo.new_segment) {
                        id++;
                        nscore += id * input.mults[y2][x2];
                    }
                    for (int i = dinfo.end; i < (int)move_seq.size(); i++) {
                        id++;
                        auto [y2, x2] = move_seq[i];
                        nscore += id * input.mults[y2][x2];
                    }
                    if (score < nscore) { // new path found!
                        dump(score, nscore, dinfo.old_segment.size(), dinfo.new_segment.size());
                        dinfo.found = true;
                        dinfo.new_score = nscore;
                    }
                }
                return;
            }

            used[y][x] = true;
            dinfo.new_segment.emplace_back(y, x);

            int dir = input.arrows[y][x];
            for (int k = 1;; k++) {
                int ny = y + dy[dir] * k, nx = x + dx[dir] * k;
                if (input.arrows[ny][nx] == -1) break;
                if (used[ny][nx]) continue;
                dfs(dinfo, ny, nx);
                if (dinfo.found) {
                    used[y][x] = false;
                    return;
                }
            }

            dinfo.new_segment.pop_back();
            used[y][x] = false;
        }

        // [begin, end) を破壊して繋ぎ直す
        DFSInfo check_mutate(int begin, int end) {
            assert(begin >= 2); // 中間の破壊のみ考慮
            assert(begin < end);
            assert(end < (int)move_seq.size());

            DFSInfo dinfo;
            dinfo.begin = begin;
            dinfo.end = end;

            for (int i = begin; i < end; i++) {
                const auto& [y, x] = move_seq[i];
                assert(used[y][x]);
                used[y][x] = false;
                dinfo.old_segment.emplace_back(y, x);
            }

            auto [sy, sx] = move_seq[begin - 1];
            std::tie(dinfo.ty, dinfo.tx) = move_seq[end];
            assert(used[dinfo.ty][dinfo.tx]);
            used[dinfo.ty][dinfo.tx] = false;
            
            // (sy, sx) -> (ty, tx) のルートを探す
            int dir = input.arrows[sy][sx];
            for (int k = 1;; k++) {
                int y = sy + dy[dir] * k, x = sx + dx[dir] * k;
                if (input.arrows[y][x] == -1) break;
                if (used[y][x]) continue;
                dfs(dinfo, y, x);
                if (dinfo.found) break;
            }

            for (int i = begin; i < end; i++) {
                const auto& [y, x] = move_seq[i];
                assert(!used[y][x]);
                used[y][x] = true;
            }
            assert(!used[dinfo.ty][dinfo.tx]);
            used[dinfo.ty][dinfo.tx] = true;

            return dinfo;
        }

        void accept_mutate(const DFSInfo& dinfo) {
            std::vector<int> new_mults;
            for (const auto& [y, x] : dinfo.new_segment) {
                new_mults.push_back(input.mults[y][x]);
            }
            for (int i = dinfo.begin; i < dinfo.end; i++) {
                const auto& [y, x] = move_seq[i];
                assert(used[y][x]);
                used[y][x] = false;
            }
            for (const auto& [y, x] : dinfo.new_segment) {
                assert(!used[y][x]);
                used[y][x] = true;
            }
            {
                auto it = move_seq.erase(move_seq.begin() + dinfo.begin, move_seq.begin() + dinfo.end);
                move_seq.insert(it, dinfo.new_segment.begin(), dinfo.new_segment.end());
            }
            { // いらんかも
                auto it = mult_seq.erase(mult_seq.begin() + dinfo.begin, mult_seq.begin() + dinfo.end);
                mult_seq.insert(it, new_mults.begin(), new_mults.end());
            }
            score = dinfo.new_score;
        }

        bool sequential_search(const int max_length, Timer& timer) {
            for (int begin = 2; begin < (int)move_seq.size() - 2; begin++) {
                for (int end = begin + 1; end < (int)move_seq.size(); end++) {
                    if (timer.elapsed_ms() > 9500) return false;
                    if (end - begin > max_length) continue;
                    auto dinfo = check_mutate(begin, end);
                    if (dinfo.found) {
                        accept_mutate(dinfo);
                        return true;
                    }
                }
            }
            return false;
        }

    };

    std::pair<int, std::vector<std::pair<int, int>>> solve(
        const Input& input,
        std::vector<std::pair<int, int>> moves,
        Timer& timer
    ) {

        State state(input, moves);
        while (state.sequential_search(10, timer)) {}

        moves = state.move_seq;
        moves.erase(moves.begin());

        return { state.score, moves };
    }

}


int main(int argc, char** argv) {

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    const bool LOCAL_MODE = argc > 1 && std::string(argv[1]) == "local";
    const int seed = 1;

    const auto input = [&]() {
        if (LOCAL_MODE) {
            return Input::load(seed);
        }
        return Input::load(std::cin);
    }();

    int best_score = 0;
    std::vector<std::pair<int, int>> best_moves;
    {
        auto [score, moves] = ReversedDFS::solve(input, 3000);
        if (chmax(best_score, score)) {
            best_moves = moves;
        }
    }
    {
        auto [score, moves] = PathMutation::solve(input, best_moves, timer);
        if (chmax(best_score, score)) {
            best_moves = moves;
        }
    }

    if (LOCAL_MODE) {
        std::ofstream out(format("../../tester/out/%d.out", seed));
        out << best_moves.size() << '\n';
        for (const auto& [y, x] : best_moves) out << y - 1 << ' ' << x - 1 << '\n';
    }
    else {
        std::ostream& out = std::cout;
        out << best_moves.size() << '\n';
        for (const auto& [y, x] : best_moves) out << y - 1 << ' ' << x - 1 << '\n';
    }

    dump(timer.elapsed_ms());

    return 0;
}