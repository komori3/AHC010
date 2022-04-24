#define _CRT_NONSTDC_NO_WARNINGS
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
//#include <atcoder/all>
//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/multiprecision/cpp_bin_float.hpp>
#ifdef _MSC_VER
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
int __builtin_clz(unsigned int n)
{
    unsigned long index;
    _BitScanReverse(&index, n);
    return 31 - index;
}
int __builtin_ctz(unsigned int n)
{
    unsigned long index;
    _BitScanForward(&index, n);
    return index;
}
namespace std {
    inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); }
}
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) {
    std::string s(v.size(), ' ');
    for (int i = 0; i < v.size(); i++) s[i] = v[i] + '0';
    os << s;
    return os;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 3.0e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 3.0e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};

/* rand */
struct Xorshift {
    uint64_t x = 88172645463325252LL;
    void set_seed(unsigned seed, int rep = 100) { x = uint64_t((seed + 1) * 10007); for (int i = 0; i < rep; i++) next_int(); }
    unsigned next_int() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    unsigned next_int(unsigned mod) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % mod; }
    unsigned next_int(unsigned l, unsigned r) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % (r - l + 1) + l; } // inclusive
    double next_double() { return double(next_int()) / UINT_MAX; }
} rnd;

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

using ll = long long;
using ld = double;
//using ld = boost::multiprecision::cpp_bin_float_quad;
using pii = std::pair<int, int>;
using pll = std::pair<ll, ll>;
using namespace std;



constexpr int N = 30;
constexpr int inf = INT_MAX / 8;
constexpr int di[] = { 0, -1, 0, 1 };
constexpr int dj[] = { -1, 0, 1, 0 };
const string d2c = "LURD";
int c2d[256];

constexpr int ROTATE[] = { 1,2,3,0,5,4,7,6 };
constexpr int TO[8][4] = {
    {1, 0, -1, -1},
    {3, -1, -1, 0},
    {-1, -1, 3, 2},
    {-1, 2, 1, -1},
    {1, 0, 3, 2},
    {3, 2, 1, 0},
    {2, -1, 0, -1},
    {-1, 3, -1, 1},
};



struct Input {
    int T[N][N];
    Input(int seed) {
        Xorshift rnd; rnd.set_seed(seed);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int r = rnd.next_int(4);
                if (r < 1) {
                    T[i][j] = rnd.next_int(4);
                }
                else if (r < 3) {
                    T[i][j] = rnd.next_int(4, 5);
                }
                else {
                    T[i][j] = rnd.next_int(6, 7);
                }
            }
        }
    }
    Input(istream& in) {
        string buf;
        for (int i = 0; i < N; i++) {
            in >> buf;
            for (int j = 0; j < N; j++) T[i][j] = buf[j] - '0';
        }
    }
    string stringify() const {
        string res;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                res += char(T[i][j] + '0');
            }
            res += '\n';
        }
        return res;
    }
};

struct Output {
    int R[N][N];
    Output() {
        memset(R, 0, sizeof(int) * N * N);
    }
    string stringify() const {
        string res;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                res += char(R[i][j] + '0');
            }
        }
        return res;
    }
    static Output generate_random(Xorshift& rnd) {
        Output out;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                out.R[i][j] = rnd.next_int(4);
            }
        }
        return out;
    }
};



double compute_score(const Input& in, const Output& out) {
    int tiles[N][N];
    memcpy(tiles, in.T, sizeof(int) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < out.R[i][j]; k++) {
                tiles[i][j] = ROTATE[tiles[i][j]];
            }
        }
    }
    vector<int> ls;
    bool used[N][N][4] = {};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int d = 0; d < 4; d++) {
                if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
                    int i2 = i, j2 = j, d2 = d, length = 0;
                    while (!used[i2][j2][d2]) {
                        if (TO[tiles[i2][j2]][d2] == -1) break;
                        length++;
                        used[i2][j2][d2] = true;
                        d2 = TO[tiles[i2][j2]][d2];
                        used[i2][j2][d2] = true;
                        i2 += di[d2];
                        j2 += dj[d2];
                        if (i2 < 0 || j2 < 0 || i2 >= N || j2 >= N) break;
                        d2 = (d2 + 2) & 3;
                    }
                    if (i == i2 && j == j2 && d == d2) {
                        ls.push_back(length);
                    }
                }
            }
        }
    }
    if (ls.size() <= 1) return 0;
    sort(ls.rbegin(), ls.rend());
    return ls[0] * ls[1];
}

struct Solver {

    Timer timer;
    Xorshift rnd;
    Input input;
    Output output;

    Solver(const Input& input) : input(input) {}

    double compute_score(const Output& out) const {
        int tiles[N][N];
        memcpy(tiles, input.T, sizeof(int) * N * N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < out.R[i][j]; k++) {
                    tiles[i][j] = ROTATE[tiles[i][j]];
                }
            }
        }
        vector<int> ls;
        bool used[N][N][4] = {};
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int d = 0; d < 4; d++) {
                    if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
                        int i2 = i, j2 = j, d2 = d, length = 0;
                        while (!used[i2][j2][d2]) {
                            if (TO[tiles[i2][j2]][d2] == -1) break;
                            length++;
                            used[i2][j2][d2] = true;
                            d2 = TO[tiles[i2][j2]][d2];
                            used[i2][j2][d2] = true;
                            i2 += di[d2];
                            j2 += dj[d2];
                            if (i2 < 0 || j2 < 0 || i2 >= N || j2 >= N) break;
                            d2 = (d2 + 2) & 3;
                        }
                        if (i == i2 && j == j2 && d == d2) {
                            ls.push_back(length);
                        }
                    }
                }
            }
        }
        if (ls.size() <= 1) return 0;
        sort(ls.rbegin(), ls.rend());
        return ls[0] * ls[1];
    }

    std::pair<double, string> solve() {

        int prev_score = compute_score(output);
        int loop = 0;
        while (timer.elapsed_ms() < 1900) {
            auto new_output = Output::generate_random(rnd);
            int score = compute_score(new_output);
            if (prev_score < score) {
                dump(score);
                prev_score = score;
                output = new_output;
            }
            loop++;
        }
        dump(loop);

        return { prev_score, output.stringify() };
    }
};


#ifdef _MSC_VER
void batch_test(int seed_begin = 0, int num_seed = 100) {

    constexpr int batch_size = 10;
    int seed_end = seed_begin + num_seed;
    
    vector<double> scores(num_seed, 0.0);
    concurrency::critical_section mtx;
    for (int batch_begin = seed_begin; batch_begin < seed_end; batch_begin += batch_size) {
        int batch_end = std::min(batch_begin + batch_size, seed_end);
        concurrency::parallel_for(batch_begin, batch_end, [&mtx, &scores](int seed) {
            std::ifstream ifs(format("tools/in/%04d.txt", seed));
            std::istream& in = ifs;
            std::ofstream ofs(format("tools/out/%04d.txt", seed));
            std::ostream& out = ofs;

            Input tc(in);

            Solver solver(tc);
            auto [score, ans] = solver.solve();

            {
                mtx.lock();
                scores[seed] = score;
                cerr << seed << ": " << score << '\n';
                mtx.unlock();
            }

            out << ans << endl;
            ifs.close();
            ofs.close();
        });
    }

    dump(std::accumulate(scores.begin(), scores.end(), 0.0));
}
#endif

int main(int argc, char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    c2d['L'] = 0; c2d['U'] = 1; c2d['R'] = 2; c2d['D'] = 3;

#if 0
    //batch_test();
    batch_test(0, 100);
#else
    std::istream& in = cin;
    std::ostream& out = cout;

    Input input(in);

    Output output;

    Solver solver(input);

    auto [score, ans] = solver.solve();

    dump(score, ans);

    out << ans << endl;

#endif

    return 0;
}