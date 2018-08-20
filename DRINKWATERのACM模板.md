# DRINKWATERのACM模板

## 数学

### 快速幂

```cpp
ll qpow(ll a, ll b) {
    ll ret = 1;
    while (b) {
        if (b & 1) ret = ret * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return ret;
}
```

### 扩展欧几里得

```cpp
void exgcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return;
    }
    exgcd(b, a % b, x, y);
    ll tmp = x;
    x = y;
    y = tmp - a / b * y;
}
```

### 欧拉函数

```cpp
ll phi(int a) {
    ll ret = a;
    for (ll i = 2; i * i <= a; i++) {
        if (a % i == 0) {
            ret = ret / i * (i - 1);
            while (a % i == 0) a /= i;
        }
    }
    if (a > 1) ret = ret / a * (a - 1);
    return ret;
}
```

### 线性筛

//prime 质数表，fac 最小质因子，phi 欧拉函数，mu 莫比乌斯函数

```cpp
void filter() {
    cnt = 0;
    memset(f, 0, sizeof(f));
    mu[1] = 1;
    phi[1] = 1;
    for (ll i = 2; i < n; i++) {
        if (!f[i]) {
            prime[cnt++] = i;
            fac[i] = i;
            phi[i] = i - 1;
            mu[i] = -1;
        }
        for (int j = 0; j < cnt; j++) {
            if (i * prime[j] >= n) break;
            f[i * prime[j]] = true;
            fac[i * prime[j]] = prime[j];
            if (i % prime[j] == 0) {
                phi[i * prime[j]] = phi[i] * prime[j];
                mu[i * prime[j]] = 0;
                break;
            } else {
                phi[i * prime[j]] = phi[i] * (prime[j] - 1);
                mu[i * prime[j]] = -mu[i];
            }
        }
    }
}
```

### 逆元

#### 费马小定理

```cpp
ll inv(ll a, ll mod) {
    return qpow(a, mod - 2);
}
```

#### 递推

```cpp
int[1] = 1;
for(int i = 2; i <= n; i++) {
    inv[i] = (mod - mod / i) * inv[mod % i] % mod;
}
```

### Lucas定理

#### 带预处理

```cpp
void init(ll mod) {
    fac[0] = 1;
    fac[1] = 1;
    facinv[0] = 1;
    facinv[1] = 1;
    for (int i = 2; i < n; i++) {
        fac[i] = fac[i - 1] * i % mod;
        facinv[i] = facinv[i - 1] * inv[i] % mod;
    }
}

ll C(ll n, ll m, ll mod) {
    return fac[n] * facinv[m] % mod * facinv[n - m] % mod;
}

ll Lucas(ll n, ll m, ll mod) {
    ll ret = 1;
    while (n || m) {
        ans = ans * C(n % mod, m % mod, mod) % mod;
        n /= mod;
        m /= mod;
    }
    return ret;
}
```

#### 不带预处理

```cpp
ll C(ll n, ll m) {
    if (m > n) return 0;
    ll ans = 1;
    for (int i = 1; i <= m; i++) {
        ll a = (n + i - m) % p;
        ll b = i % p;
        ans = ans * (a * qpow(b, p - 2) % p) % p;
    }
    return ans;
}

ll Lucas(ll n, ll m) {
    if (m == 0) return 1;
    return C(n % p, m % p) * Lucas(n / p, m / p) % p;
}
```

### 组合数单行递推

$C_{n,0} = 1$

$C_{n,i} = C_{n,i - 1} * \frac{n - i + 1}{i}$

### 中国剩余定理

```cpp
int CRT(int a[], int m[], int n) {
    int M = 1, ret = 0;
    for (int i = 1; i <= n; i++) M *= m[i];
    for (int i = 1; i <= n; i++) {
        int x, y;
        int Mi = M / m[i];
        exgcd(Mi, m[i], x, y);
        ret = (ret + Mi * x * a[i]) % M;
    }
    if (ret < 0) ret += M;
    return ret;
}
```

### 高斯消元

```cpp
#define eps 1e-9
const int MAXN = 220;
double a[MAXN][MAXN], x[MAXN]; // 方程的左边的矩阵和等式右边的值，求解之后 x 存的就是结果
int equ, var; // 方程数和未知数个数
/*
 * 返回 0 表示无解，1 表示有解
 */
int Gauss() {
    int i, j, k, col, max_r;
    for (k = 0, col = 0; k < equ && col < var; k++; col++) {
        max_r = k;
        for (i = k + 1; i < equ; i++)
            if (fabs(a[i][col]) > fabs(a[max_r][col]))
                max_r = i;
        if (fabs(a[max_r][col]) < eps) return 0;
        if (k != max_r) {
            for (j = col; j < var; j++)
                swap(a[k][j], a[max_r][j]);
            swap(x[k], x[max_r]);
        }
        x[k] /= a[k][col];
        for (j = col + 1; j < var; j++)
            a[k][j] /= a[k][col];
        a[k][col] = 1;
        for (i = 0; i < equ; i++)
            if (i != k) {
                x[i] -= x[k] * a[i][k];
                for (j = col + 1; j < var; j++)
                    a[i][j] -= a[k][j] * a[i][col];
                a[i][col] = 0;
            }
    }
    return 1;
}
```

### 质数计数

| $n$ | $n$ 以内质数个数 |
| - | - |
| $10^{0}$ | 0 |
| $10^{1}$ | 4 |
| $10^{2}$ | 25 |
| $10^{3}$ | 168 |
| $10^{4}$ | 1229 |
| $10^{5}$ | 9592 |
| $10^{6}$ | 78498 |
| $10^{7}$ | 664579 |
| $10^{8}$ | 5761455 |
| $10^{9}$ | 50847534 |
| $10^{10}$ | 455052511 |
| $10^{11}$ | 4118054813 |
| $10^{12}$ | 37607912018 |
| $10^{13}$ | 346065536839 |
| $10^{14}$ | 3204941750802 |
| $10^{15}$ | 29844570422669 |
| $10^{16}$ | 279238341033925 |
| $10^{17}$ | 2623557157654233 |
| $10^{18}$ | 24739954287740860 |

### FFT

```cpp
struct plex {
    double x, y;
    plex (double _x = 0.0, double _y = 0.0) : x (_x), y (_y) {}
    plex operator +(const plex &a) const {
        return plex(x + a.x, y + a.y);
    }
    plex operator -(const plex &a) const {
        return plex(x - a.x, y - a.y);
    }
    plex operator *(const plex &a) const {
        return plex(x * a.x - y * a.y, x * a.y + y * a.x);
    }
};

void change(plex *y, int len) {
    int i, j, k;
    for (i = 1, j = len / 2; i < len - 1; i++) {
        if (i < j) swap(y[i], y[j]);
        k = len / 2;
        while (j >= k) {
            j -= k;
            k /= 2;
        }
        if (j < k) j += k;
    }
}

void fft(plex y[], int len, int on) {
    change(y, len);
    for (int h = 2; h <= len; h <<= 1) {
        plex wn(cos(-on * 2 * pi / h), sin(-on * 2 * pi / h));
        for (int j = 0; j < len; j += h) {
            plex w(1, 0);
            for (int k = j; k < j + h / 2; k++) {
                plex u = y[k];
                plex t = w * y[k + h / 2];
                y[k] = u + t;
                y[k + h / 2] = u - t;
                w = w * wn;
            }
        }
    }
    if (on == -1)
        for (int i = 0; i < len; i++)
            y[i].x /= len;
}
```

### NTT

```cpp
ll x1[MAXN], x2[MAXN];
void change(ll y[], int len) {
    for (int i = 1, j = len / 2; i < len - 1; i++) {
        if (i < j) swap(y[i], y[j]);
        int k = len / 2;
        while (j >= k) {
            j -= k;
            k /= 2;
        }
        if (j < k) j += k;
    }
}

void ntt(ll y[], int len, int on) {
    change(y, len);
    for (int h = 2; h <= len; h <<= 1) {
        ll wn = qpow(G, (mod - 1) / h);
        if (on == -1) wn = qpow(wn, mod - 2);
        for (int j = 0; j < len; j += h) {
            ll w = 1;
            for (int k = j; k < j + h / 2; k++) {
                ll u = y[k];
                ll t = w * y[k + h / 2] % mod;
                y[k] = (u + t) % mod;
                y[k + h / 2] = (u - t + mod) % mod;
                w = w * wn % mod;
            }
        }
    }
    if (on == -1) {
        ll t = qpow(len, mod - 2);
        for (int i = 0; i < len; i++)
            y[i] = y[i] * t % mod;
    }
}
```

### 数列

#### Catalan 数列

A000108

1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012, 742900, 2674440, 9694845, 35357670, 129644790, 477638700, 1767263190, 6564120420, 24466267020, 91482563640, 343059613650, 1289904147324, 4861946401452, 18367353072152, 69533550916004, 263747951750360, 1002242216651368, 3814986502092304

$a_{n} = \frac{1}{n + 1} {2n \choose n}$

$a_{n} = {2n \choose n} - {2n \choose n + 1}$

$a_{i} = \frac{i * 4 - 2}{i + 1}a_{i - 1}$ ($a_{1} = 1$)

#### 九连环数列

A000975
0, 1, 2, 5, 10, 21, 42, 85, 170, 341, 682, 1365, 2730, 5461, 10922, 21845, 43690, 87381, 174762, 349525, 699050, 1398101, 2796202, 5592405, 11184810, 22369621, 44739242, 89478485, 178956970, 357913941, 715827882, 1431655765, 2863311530

$a_{2n} = 2 * a_{2n - 1}$

$a_{2n + 1} = 2 * a_{2n} + 1$

${a_{n}}$ 是第n个二进制表示中没有连续相同数字的数

#### 第二类 Stirling 数

$a_{i, j} = a_{i - 1, j - 1} + j * a_{i - 1, j}$

#### 接近正三角形的海伦三角形

A003500

2, 4, 14, 52, 194, 724, 2702, 10084, 37634, 140452, 524174, 1956244, 7300802, 27246964, 101687054, 379501252, 1416317954, 5285770564, 19726764302, 736212866644, 274758382274, 1025412242452, 3826890587534, 14282150107684

### 异或线性基

```cpp
void cal() {
    for (int i = 0; i < n; i++)
        for (int j = MAX_BASE; j >= 0; j--)
            if (a[i] >> j & 1) {
                if (b[j]) a[i] ^= b[j];
                else {
                    b[j] = a[i];
                    for (int k = j - 1; k >= 0; k--)
                        if (b[k] && (b[j] >> k & 1)) b[j] ^= b[k];
                    for (int k = j + 1; k <= MAX_BASE; k++)
                        if (b[k] >> j & 1) b[k] ^= b[j];
                    break;
                }
            }
}
```

## 几何

### 海伦公式

$p = \frac{a + b + c}{2}$

$S = \sqrt{p(p - a)(p - b)(p - c)}$

## 数据结构

### 树状数组

```cpp
int lowbit(int a) {
    return a & (-a);
}

void add(int x, int v) {
    while (x <= n) {
        BIT[x] += v;
        x += lowbit(x);
    }
}

int sum(int x) {
    int ret = 0;
    while (x > 0) {
        ret += BIT[x];
        x -= lowbit(x);
    }
    return ret;
}
```

### 线段树

```cpp
#define lchild x << 1, l, mid
#define rchild x << 1 | 1, mid + 1, r

void build(int x = 1, int l = 1, int r = n) {
    if (l == r) {
        st[x] = a[l];
        return;
    }
    int mid = (l + r) >> 1;
    build(lchild);
    build(rchild);
    push_up(x);
}

// 区间求和
void push_up(int x) {
    st[x] = st[x << 1] + st[x << 1 | 1];
}

void push_down(int x, int len) {
    st[x << 1] += mark[x] * (len - (len >> 1));
    mark[x << 1] += mark[x];
    st[x << 1 | 1] += mark[x] * (len >> 1);
    mark[x << 1 | 1] += mark[x];
    mark[x] = 0;
}

void query(int L, int R, int x = 1, int l = 1, int r = n) {
    if (L <= l && r <= R) return st[x];
    if (mark[x]) push_down(x, r - l + 1);
    int mid = (l + r) >> 1, ret = 0;
    if (L <= mid) ret += query(L, R, lchild);
    if (R > mid) ret += query(L, R, rchild);
    return ret;
}

void update(int L, int R, int val, int x = 1, int l = 1, int r = n) {
    if (L <= l && r <= R) {
        st[x] += val * (r - l + 1);
        mark[x] += val;
        return;
    }
    if (mark[x]) push_down(x, r - l + 1);
    int mid = (l + r) >> 1;
    if (L <= mid) update(L, R, val, lchild);
    if (R > mid) update(L, R, val, rchild);
    push_up(x);
}

// 区间求最大值
void push_up(int x) {
    st[x] = max(st[x << 1], st[x << 1 | 1]);
}

void push_down(int x) {
    st[x << 1] += mark[x];
    mark[x << 1] += mark[x];
    st[x << 1 | 1] += mark[x];
    mark[x << 1 | 1] += mark[x];
    mark[x] = 0;
}

void query(int L, int R, int x = 1, int l = 1, int r = n) {
    if (L <= l && r <= R) return st[x];
    if (mark[x]) push_down(x);
    int mid = (l + r) >> 1, ret = 0;
    if (L <= mid) ret = max(ret, query(L, R, lchild));
    if (R > mid) ret = max(ret, query(L, R, rchild));
    return ret;
}

void update(int L, int R, int val, int x = 1, int l = 1, int r = n) {
    if (L <= l && r <= R) {
        st[x] += val;
        mark[x] += val;
        return;
    }
    if (mark[x]) push_down(x);
    int mid = (l + r) >> 1;
    if (L <= mid) update(L, R, val, lchild);
    if (R > mid) update(L, R, val, rchild);
    push_up(x);
}
```

### RMQ

```cpp
void rmq_init() {
    for (int i = 0; i < n; i++) f[i][0] = a[i];
    for (int j = 1; (1 << j) <= n; j++) {
        for (int i = 0; i + (1 << j) - 1 < n; i++) {
            f[i][j] = max(f[i][j - 1], f[i + (1 << (j - 1))][j - 1]);
        }
    }
}

int rmq(int l, int r) {
    int i = 0;
    while ((1 << (i + 1)) <= r - l + 1) i++;
    return max(f[l][i], f[r - (1 << i) + 1][i]);
}
```

### 主席树

```cpp
#define maxn 111111

struct node {
    int l, r;
    int num;
} tree[maxn << 6];
int T[maxn];
int cnt;
int n, q, a[maxn];
vector<int> num;
int p[maxn];
map<int, int> gg;

void init() {
    sort(num.begin(), num.end());
    int cnt = 0;
    for (int i = 0; i < num.size(); i++) {
        if (!i || num[i] != num[i - 1]) {
            gg[num[i]] = ++cnt;
            p[cnt] = num[i];
        }
    }
}

int build_tree(int l, int r) {
    int root = cnt++;
    tree[root].num = 0;
    if (l == r) return root;
    int mid = (l + r) >> 1;
    tree[root].l = build_tree(l, mid);
    tree[root].r = build_tree(mid + 1, r);
    return root;
}

int update(int root, int pos, int val) {
    int new_root = cnt++;
    int tmp = new_root;
    tree[new_root].nu = tree[root].num + val;
    int l = 1, r = n;
    while (l < r) {
        int mid = (l + r) >> 1;
        if (pos <= mid) {
            tree[new_root].l = cnt++;
            tree[new_root].r = tree[root].r;
            new_root = tree[new_root].l;
            root = tree[root].l;
            r = mid;
        } else {
            tree[new_root].l = tree[root].l;
            tree[new_root].r = cnt++;
            new_root = tree[new_root].r;
            root = tree[root].r;
            l = mid + 1;
        }
        tree[new_root].num = tree[root].num + val;
    }
    return tmp;
}

int query(int l_root, int r_root, int k) {
    int l = 1, r = n;
    while (l < r) {
        int mid = (l + r) >> 1;
        if (tree[tree[l_root].l].num - tree[tree[r_root].l].num >= k) {
            r = mid;
            l_root = tree[l_root].l;
            r_root = tree[r_root].l;
        } else {
            l = mid + 1;
            k -= (tree[tree[l_root].l].num - tree[tree[r_root].l].num);
            l_root = tree[l_root].r;
            r_root = tree[r_root].r;
        }
    }
    return l;
}

int main() {
    int t;
    scanf("%d", &t);
    while(t--) {
        scanf("%d%d", &n, &q);
        num.clear();
        gg.clear();
        for (int i = 1; i <= m; i++) {
            scanf("%d", &a[i]);
            num.push_back(a[i]);
        }
        init();
        cnt = 0;
        T[m + 1] = build_tree(1, n);
        for (int i = n; i >= 1; i--) {
            int pos = gg[a[i]];
            T[i] = update(T[i + 1], pos, 1);
        }
        while (q--) {
            int l, r, k;
            scanf("%d%d%d", &l, &r, &k);
            printf("%d\n", p[query(T[l], T[r + 1], k)]);
        }
    }
    return 0;
}
```

## 字符串

### KMP

```cpp
void get_fail(char *p) {
    int i = -1, j = 0, lenp = strlen(p);
    f[0] = -1;
    while (j < lenp) {
        if (i < 0 || p[i] == p[j]) {
            i++;
            j++;
            f[j] = i;
        } else i = f[i];
    }
}

int kmp(char *t, char *p) {
    int i = 0, j = 0, lent = strlen(t), lenp = strlen(p), ret = 0;
    while (i < lent && j < lenp) {
        if (j < 0 || t[i] == p[j]) {
            if(j == lenp - 1) {
                ret++;
                j = f[j];
                continue;
            }
            i++;
            j++;
        } else j = f[j];
    }
    return ret;
}
```

### Manacher

```cpp
int manacher(char *p) {
    char s[MAXN * 2];
    int lenp = strlen(p), k = 0;
    s[k++] = '@';
    s[k++] = '#';
    for (int i = 0; i < lenp; i++) {
        s[k++] = p[i];
        s[k++] = '#';
    }
    s[k++] = '~';
    s[k] = 0;
    int max = 0, pos = 0, ret = 0;
    for (int i = 1; i < k; i++) {
        if (max > i) {
            len[i] = min(len[2 * pos - i], max - i);
        } else len[i] = 1;
        while(s[i + len[i]] == s[i - len[i]]) len[i]++;
        ret = max(ret, len[i]);
        if (len[i] + i > max) {
            max = len[i] + i;
            pos = i;
        }
    }
    return ret - 1;
}
```

### 字典树

```cpp
struct trie {
    int next[maxn][26], end[maxn];
    init root, cnt;

    int new_node() {
        memset(next[cnt], -1, sizeof(next[cnt]);
        end[cnt++] = 0;
        return cnt - 1;
    }

    void init() {
        cnt = 0;
        root = new_node();
    }

    void insert(char *buf) {
        int len = strlen(buf);
        int now = root;
        for (int i = 0; i < len; i++) {
            int id = buf[i] - 'a';
            if (next[now][id] == -1) {
                next[now][id] = new_node();
            }
            now = next[now][id];
        }
        end[now]++;
    }
};
```

### AC 自动机

```cpp
struct trie {
    int next[maxn][26], fail[maxn], end[maxn];
    int root, cnt;

    int new_node() {
        memset(next[cnt], -1, sizeof(next[cnt]);
        end[cnt++] = 0;
        return cnt - 1;
    }

    void init() {
        cnt = 0;
        root = new_node();
    }

    void insert(char *buf) {
        int len = strlen(buf);
        int now = root;
        for (int i = 0; i < len; i++) {
            int id = buf[i] - 'a';
            if (next[now][id] == -1) {
                next[now][id] = new_node();
            }
            now = next[now][id];
        }
        end[now]++;
    }

    void build() {
        queue<int> q;
        fail[root] = root;
        for (int i = 0; i < 26; i++) {
            if (next[root][i] == -1) {
                next[root[i] = root;
            } else {
                fail[next[root][i]] = root;
                q.push(next[root][i]);
            }
        }
        while (!q.empty()) {
            int now = q.front(); q.pop();
            for (int i = 0; i < 26; i++) {
                if (next[now][i] == -1) {
                    next[now][i] = next[fail[now]][i];
                } else {
                    fail[next[now][i]] = next[fail[now]][i];
                    q.push(next[now][i]);
                }
            }
        }
    }

    int query(char *buf) {
        int len = strlen(buf);
        int now = root;
        int res = 0;
        for (int i = 0; i < len; i++) {
            int id = buf[i] - 'a';
            now = next[now][id];
            int tmp = now;
            while (tmp != root) {
                if (end[tmp]) {
                    res += end[tmp];
                    end[tmp] = 0;
                }
                tmp = fail[tmp];
            }
        }
        return res;
    }
}ac;
```

### 后缀数组

```cpp
// m 表示字符中最大的值
int t1[maxn], t2[maxn], c[maxn];

bool cmp(int *r, int a, int b, int l) {
    return r[a] == r[b] && r[a + l] == r[b + l];
}

void da(int str[], int sa[], int rank[], int height[], int n, int m) {
    n++;
    int i, j, p, *x = t1, *y = t2;
    // 第一轮基数排序，如果 s 的最大值很大，可改为快速排序
    for (i = 0; i < m; i++) c[i] = 0;
    for (i = 0; i < n; i++) c[x[i] = str[i]]++;
    for (i = 1; i < m; i++) c[i] += c[i - 1];
    for (i = n - 1; i >= 0; i--) sa[--c[x[i]]] = i;
    for (j = 1; j <= n; j <<= 1) {
        p = 0;
        // 直接利用 sa 数组排序第二关键字
        for (i = n - j; i < n; i++) y[p++] = i; // 后面的 j 个数第二关键字为空的最小
        for (i = 0; i < n; i++) if (sa[i] >= j) y[p++] = sa[i] - j;
        // 这样数组 y 保存的就是按照第二关键字排序的结果
        // 基数排序第一关键字
        for (i = 0; i < m; i++) c[i] = 0;
        for (i = 0; i < n; i++) c[x[y[i]]]++;
        for (i = 1; i < m; i++) c[i] += c[i - 1];
        for (i = n - 1; i >= 0; i--) sa[--c[x[y[i]]] = y[i];
        // 根据 sa 和 x 数组计算新的 x 数组
        swap(x, y);
        p = 1;
        x[sa[0]] = 0;
        for (i = 1; i < n; i++)
            x[sa[i]] = cmp(y, sa[i - 1], sa[i], j) ? p - 1 : p++;
        if (p >= n) break;
        m = p; // 下次基数排序的最大值
    }
    int k = 0;
    n--;
    for (i = 0; i <= n; i++) rank[sa[i]] = i;
    for (i = 0; i < n; i++) {
        if (k) k--;
        j = sa[rank[i] - 1];
        while (str[i + k] == str[j + k]) k++;
        height[rank[i]] = k;
    }
}
```

## 图论

__双向边开两倍！！！双向边开两倍！！！双向边开两倍！！！__

### 链式前向星

```cpp
struct Edge {
    int next, v, w;
};

void add(int u, int v, int w) {
    edge[cnt].w = w;
    edge[cnt].v = v;
    edge[cnt].next = head[u];
    head[u] = cnt++;
}
// head[] 数组初始化为 -1
```

### 最短路

#### 单源最短路

##### Dijsktra

###### 朴素写法

```cpp
for (int i = 1; i <= n; i++) dis[i] = INF;
dis[s] = 0;
for (int i = 1; i <= n; i++) {
    int x, mind = INF;
    for (int j = 1; j <= n; j++)
        if (!f[j] && dis[j] < mind) {
            mind = dis[j];
            x = j;
        }
    f[x] = true;
    for (int j = 1; j <= n; j++) dis[j] = min(dis[j], dis[x] + w[x][j]);
}
```

###### 堆优化

```cpp
bool vis[maxn];
ll d[maxn];

struct E {
    int v; ll w;
    bool operator <(const E &_) const {
        return w > _.w;
    }
};

void dis(int s, int n) {
    for (int i = 1; i <= n; i++) {
        d[i] = INF;
        vis[i] = 0;
    }
    vis[s] = 1; d[s] = 0;
    priority_queue<E> q;
    for (int i = head[s]; i + 1; i = edge[i].next) {
        q.push((E) {edge[i].v, edge[i].w});
    }
    while (!q.empty()) {
        E tmp = q.top(); q.pop();
        int v = tmp.v; ll w = tmp.w;
        if (vis[tmp.v]) continue;
        d[v] = w;
        vis[v] = 1;
        for (int i = head[v]; i + 1; i = edge[i].next) {
            if (!vis[edge[i].v])
                q.push((E) {edge[i].v, d[v] + edge[i].w});
        }
    }
}
```

##### SPFA

```cpp
// 返回 1 表示没有负环，否则表示存在负环。
// 如果要判断全图是否有负环，初始化时将全部点入队列，否则初始化源点入队。
bool vis[maxn];
int top, num[maxn];
double d[maxn];
bool spfa(int n) {
    queue<int> q;
    while (!q.empty()) q.pop();
    memset(num, 0, sizeof(num));
    for (int i = 0; i < n; i++) {
        vis[i] = 1;
        d[i] = 0;
        q.push(i);
        num[i] = 1;
    }
    while(!q.empty()) {
        int u = q.front(); q.pop();
        vis[u] = 0;
        for (int i = head[u]; i != -1; i = edge[i].next) {
            int v = edge[i].v;
            if (d[v] > d[u] + edge[i].w) {
                d[v] = d[u] + edge[i].w;
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                    if (++num[v] > n) return 0;
                }
            }
        }
    }
    return 1;
}
```

#### 多源最短路

##### Floyd

```cpp
// 初始化：mp[i][i] = 初始值 （0 或者 INF， 根据题意来）
// 如果 <i, j> 没有边，dp[i][j] = INF（小心溢出）
void floyd() {
    for (int k = 1; k <= n; k++) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                mp[i][j] = min(mp[i][j], mp[i][k] + mp[k][j]);
            }
        }
    }
}
```

### 次短路

```cpp
#include<cstdio>
#include<cstring>
#include<queue>
#include<algorithm>
#define MAXN (5000 + 10)
#define INF (5000 * 5000 * 2)
using namespace std;

typedef pair<int, int> P;

struct edge{
    int to, cost;
    edge(int tv = 0, int tc = 0):
        to(tv), cost(tc){}
};

int N, R;
vector<edge> graph[MAXN];
int dist[MAXN]; // 最短距离
int dist2[MAXN]; // 次短距离

void solve() {
    fill(dist, dist + N, INF);
    fill(dist2, dist2 + N, INF);
    // 从小到大的优先队列
    // 使用 pair 而不用 edge 结构体
    // 是因为这样我们不需要重载运算符
    // pair 是以 first 为主要关键字进行排序
    priority_queue<P, vector<P>, greater<P>> Q;
    // 初始化源点信息
    dist[0] = 0;
    Q.push(P(0, 0));
    // 同时求解最短路和次短路
    while (!Q.empty()) {
        P p = Q.top(); Q.pop();
        // first 为 s -> to 的距离，second 为 edge 结构体的 to
        int v = p.second, d = p.first
        // 当取出的值不是当前最短距离或次短距离，就舍弃他
        if (dist2[v] < d) continue;
        for (unsigned i = 0; i < graph[v].size(); i++) {
            edge &e = graph[v][i];
            int d2 = d + e.cost;
            if (dist[e.to] > d2) {
                swap(dist[e.to], d2);
                Q.push(P(dist[e.to], e.to));
            }
            if (dist2[e.to] > d2 && dist[v] < d2) {
                dist2[e.to] = d2;
                Q.push(P(dist2[e.to], e.to));
            }
        }
    }
    printf("%d\n", dist2[N - 1]);
}

int main() {
    int A, B, D;
    scanf("%d%d", &N, &R);
    for (int i = 0; i < R; i++) {
        scanf("%d%d%d", &A, &B, &D);
        graph[A - 1].push_back(edge(B - 1, D));
        graph[B - 1].push_back(edge(A - 1, D));
    }
    solve();
    return 0;
}
```

### 最小生成树

#### Prim

```cpp
const int INF = 0x3f3f3f3f;
const int MAXN = 110;
bool vis[MAXN];
int lowc[MAXN];
int Prim(int cost[][MAXN], int n) { // 点是 0 ~ n - 1
    int ans = 0;
    memset(vis, false, sizeof(vis));
    vis[0] = true;
    for (int i = 1; i < n; i++) lowc[i] = cost[0][i];
    for (int i = 1; i < n; i++) {
        int minc = INF;
        int p = -1;
        for (int j = 0; j < n; j++)
            if (!vis[j] && minc > lowc[j]) {
                minc = lowc[j];
                p = j;
            }
            if (minc == INF) return -1;
            ans += minc;
            vis[p] = true;
            for (int j = 0; j < n; j++)
                if (!vis[j] && lowc[j] > cost[p][j])
                    lowc[j] = cost[p][j];
    }
    return ans;
}
```

#### Kruskal

```cpp
const int MAXN = 110; // 最大点数
const int MAXM = 10000; // 最大边数
int F[MAXN]; // 并查集使用
struct Edge {
    int u, v, w;
} edge[MAXM]; // 存储边的信息，包括起点/终点/权值
int tol; // 边数，加边前赋值为 0
void addedge(int u, int v, int w) {
    edge[tol].u = u;
    edge[tol].v = v;
    edge[tol++].w = w;
}
bool cmp(Edge a, Edge b) { // 排序函数，将边按照权值从小到大排序
    return a.w < b.w;
}
int find(int x) {
    if (F[x] == -1) return x;
    else return F[x] = find(F[x]);
}
int Kruskal(int n) { // 传入点数，返回最小生成树的权值，如果不连通返回 -1
    memset(F, -1, sizeof(F));
    sort(edge, edge + tol, cmp);
    int cnt = 0; // 计算加入的边数
    int ans = 0;
    for (int i = 0; i < tol; i++) {
        int u = edge[i].u;
        int v = edge[i].v;
        int w = edge[i].w;
        int t1 = find(u);
        int t2 = find(v);
        if (t1 != t2) {
            ans += w;
            F[t1] = t2;
            cnt++;
        }
        if (cnt == n - 1) break;
    }
    if (cnt < n - 1) return -1; // 不连通
    else return ans;
}
```

### 最近公共祖先

#### 在线 RMQ 查询

```cpp
void dfs(int k, int depth, int par) {
    D[k] = depth;
    m++;
    E[m] = k;
    R[k] = m;
    for (vector<int>::iterator it = edge[k].begin(); it != edge[k].end(); it++) {
        if (*it == par) continue;
        dfs(*it, depth + 1, k);
        m++;
        E[m] = k;
    }
}

void init_RMQ() {
    g[0] = -1;
    for (int i = 1; i <= m; i++) {
        g[i] = ((i & (i - 1)) == 0) ? g[i - 1] + 1 : g[i - 1];
        f[i][0] = i;
    }
    for (int j = 1; j <= g[m]; j++)
        for (int i = 1; i + (1 << j) - 1 <= m; i++)
            dp[i][j] = D[E[f[i][j - 1]]] < D[E[f[i + (1 << (j - 1))][j - 1]]] ? f[i][j - 1] : f[i + (1 << (j - 1))][j - 1];
}

int RMQ(int x, int y) {
    if(x > y) swap(x, y);
    int k = g[y - x + 1];
    return D[E[f[x][k]]] <= D[E[f[y - (1 << k) + 1][k]]] ? f[x][k] : f[y - (1 << k) + 1][k];
}

int LCA(int u, int v) {
    return E[RMQ(R[u], R[v])];
}
```

### 割点和桥

```cpp
struct node {
    int u, v, next;
    bool is; // 这条边是不是桥
} edge[maxm];

int head[maxn], cnt;
int n, m;
int dfs_clock, low[maxn], pre[maxn];

void add_edge(int u, int v) {
    edge[cnt].is = 0;
    edge[cnt].u = u, edge[cnt].v = v, edge[cnt].next = head[u], head[u] = cnt++;
}

void find_cnt(int u, int father) {
    pre[u] = low[u] = ++dfs_clock;
    for (int i = head[u]; i + 1; i = edge[i].next) {
        int v = edge[i].v;
        if (v == father) continue;
        if (!pre[v]) {
            find_cut(v, u);
            if (low[v] <= pre[u]) {}
            else {
                edge[i].is = edge[i ^ 1].is = 1;
            }
            low[u] = min(low[u], low[v]);
        } else {
            low[u] = min(low[u], pre[v]);
        }
    }
}

int main() {
    ios::sync_with_stdio(0);
    cin >> n >> m;
    memset(head, -1, sizeof(head));
    cnt = 0;
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        add_edge(u, v);
        add_edge(v, u);
    }
    dfs_clock = 0;
    for (int i = 1; i <= n; i++)
        if (!pre[i]) find_cut(i, 0);
    return 0;
}
```

### 强连通分量

```cpp
int head[2000], edg_cnt = 0, dfn[N], low[N], sccno[N], step, scc_cnt = 0, ip;
// dfn[] 表示深搜的步数，low[u] 表示 u 或 u 的子树能够追溯到的最早的栈中节点的次序号
// sccno[] 缩点数组，表示某个点对应的缩点值
// scc_cnt 强连通分量个数
vector<int> scc[N]; // 得出来的缩点，scc[i] 里面存 i 这个缩点具体缩了哪些点
stack<int> s;

void dfs(int u) {
    dfn[u] = low[u] = ++step;
    s.push(u);
    for (int k = head[u]; k != -1; k =edge[k].next) {
        int v = edge[k].v;
        if (!dfn[v]) {
            dfs(v);
            low[u] = min(low[u], low[v]);
        } else if (!sccno[v]) {
            low[u] = min(low[u], low[v]);
        }
    }
    if (low[u] == dfn[u]) {
        scc_cnt += 1;
        scc[scc_cnt].clear();
        while (true) {
            int x = s.top();
            s.pop();
            if (sccno[x] != scc_cnt) scc[scc_cnt].push_back(x);
            sccno[x] = scc_cnt;
            if (x == u) break;
        }
    }
}

void tarjan(int n) {
    memset(sccno, 0, sizeof(sccno));
    memset(dfn, 0, sizeof(dfn));
    step = scc_cnt = 0;
    for (int i = 1; i <= n; i++) {
        if (!dfn[i]) dfs(i);
    }
}
```

### 网络流

#### 最大流

```cpp
// 用 #define type xxx 替代流的类型 (int, double, ll)
// 极其关键：边数开两倍
int n;
int s, t;

struct Edge {
    int from, to, next;
    type cap, flow;
    void get(int u, int a, int b, type c, type d) {
        from = u; to = a; next = b; cap = c; flow = d;
    }
} edge[maxm];

int tol;
int head[maxn];
int gap[maxn], dep[maxn], pre[maxn], cur[maxn];

void init() {
    tol = 0;
    memset(head, -1, sizeof(head));
}

void add_edge(int u, int v, type w, type rw = 0) {
    edge[tol].get(u, v, head[u], w, 0); head[u] = tol++;
    edge[tol].get(v, u, head[v], rw, 0); head[v] = tol++;
}

type sap(int start, int end, int N) {
    memset(gap, 0, sizeof(gap));
    memset(dep, 0, sizeof(dep));
    memcpy(cur, head, sizeof(head));
    int u = start;
    pre[u] = -1;
    gap[0] = N;
    type ans = 0;
    while (dep[start] < N) {
        if (u == end) {
            type Min = INF;
            for (int i = pre[u]; i != -1; i = pre[edge[i ^ 1].to])
                if (Min > edge[i].cap - edge[i].flow)
                    Min = edge[i].cap - edge[i].flow;
            for (int i = pre[u]; i != -1; i = pre[edge[i ^ 1].to]) {
                edge[i].flow += Min;
                edge[i ^ 1].flow -= Min;
            }
            u = start;
            ans += Min;
            continue;
        }
        bool flag = false;
        int v;
        for (int i = cur[u]; i != -1; i = edge[i].next) {
            v = edge[i].to;
            if (edge[i].cap - edge[i].flow && dep[v] + 1 == dep[u]) {
                flag = true;
                cur[u] = pre[v] = i;
                break;
            }
        }
        if (flag) {
            u = v;
            continue;
        }
        int Min = N;
        for (int i = head[u]; i != -1; i = edge[i].next)
            if (edge[i].cap - edge[i].flow && dep[edge[i].to] < Min) {
                Min = dep[edge[i].to];
                cur[u] = i;
            }
        gap[dep[u]]--;
        if (!gap[dep[u]]) return ans;
        dep[u] = Min + 1;
        gap[dep[u]]++;
        if (u != start) u = edge[pre[u] ^ 1].to;
    }
    return ans;
}
```

## 二分答案

```cpp
int find(int l, int r) {
    while (r - l > 1) {
        int mid = (l + r) >> 1;
        if (check(mid)) r = mid;
        else l = mid;
    }
    int ret = (check(l) ? l : r);
    if (!check(r)) return -1; else return ret;
}
```

## 问题实例

### 单调栈求子 1 矩阵数

```cpp
for (int i = 1; i <= n; i++) {
    head = 0;
    sum = 0;
    for (int j = 1; j<= m; j++) {
        if (f[i][j]) h[j]++; else h[j] = 0;
        tmp = 1;
        while (head && val[head] >= h[j]) {
            tmp += cnt[head];
            sum -= val[head] * cnt[head];
            head--;
        }
        head++;
        val[head] = h[j];
        cnt[head] = tmp;
        sum += val[head] * cnt[head];
    }
}
```

### 最小割的最少割边数

建边的时候每条边权 w = w * (E + 1) + 1
这样得到最大流 maxflow / (E + 1)，最少割边数 maxflow % (E + 1)

### 求柱形图的最大面积

```cpp
#include<iostream>
#include<cstdio>
#include<algorithm>
using namespace std;

typedef long long ll;
const int MAXN = 1e5 + 10;
const int INF = 0x3f3f3f3f;

ll h[MAXN];
ll l[MAXN], r[MAXN];

int main() {
    int n;
    while (scanf("%d", &n) == 1 && n) {
        h[0] = h[n + 1] = -INF;
        for (int i = 1; i <= n; i++) {
            scanf("%lld", &h[i]);
        }
        for (int i = 1; i <= n; i++) {
            l[i] = r[i] = i;
        }
        for (int i = 1; i <= n; i++) {
            while (h[l[i] - 1] >= h[i])
                l[i] = l[l[i] - 1];
        }
        for (int i = n; i >= 1; i--) {
            while (h[r[i] + 1] >= h[i])
                r[i] = r[r[i] + 1];
        }
        ll ans = 0;
        for (int i = 1; i <= n; i++) {
            ans = max(ans, (r[i] - l[i] + 1) * h[i]);
        }
        printf("%lld\n", ans);
    }
}
```

### 多重背包的 01 写法

```cpp
void completePack(int dp[], int value, int weight, int total) {
    int i;
    for (i = weight; i <= total; i++) {
        dp[i] = max(dp[i], dp[i - weight] + value);
    }
}

void ZeroOnePack(int dp[], int value, int weight, int total) {
    int i;
    for (i = total; i >= weight; i--) {
        dp[i] = max(dp[i], dp[i - weight] + value);
    }
}

// 多重背包问题 优化 一维数组 二进制思想 时间复杂度为 O(V*Σlog(n[i])
void multiPack(int dp[], int value, int weight, int amount, int total) {
    if (weight * amount > total) {
        completePack(dp, value, weight, total);
    } else {
        int k = 1;
        while (amount - k >= 0) {
            ZeroOnePack(dp, k * value, k * weight, total);
            amount -= k;
            k *= 2;
        }
        ZeroOnePack(dp, amount * value, amount * weight, total);
    }
}

int main() {
    int n, w;
    cin >> n >> w;
    int i;
    int wi, vi, ci;
    for (i = 0; i < n; i++) {
        cin >> wi >> vi >> ci;
        multiPack(dp_1, vi, wi, ci, w);
    }
    cout << dp_1[w] << endl;
    return 0;
}
```

## 高精度

### Java

#### a + b

```java
import java.math.BigInteger;
import java.util.Scanner;

public class Main {
    void solve() {
        BigInteger a, b, c;
        Scanner cin = new Scanner(System.in);
        int t = cin.nextInt();
        for (int i = 1; i <= t; i++) {
            System.out.println("Case " + i + ":");
            a = cin.nextBigInteger();
            b = cin.nextBigInteger();
            System.out.println(a + " + " + b + " = " + a.add(b));
            if (i != t) System.out.println();
        }
    }
    public static void main(String[] args) {
        Main work = new Main();
        work.solve();
    }
}
```

#### 阶乘

```java
import java.math.BigInteger;
import java.util.Scanner;

public class Main {
    int maxn = 10005;
    void solve() {
        Scanner cin = new Scanner(System.in);
        int n;
        while (cin.hasNext()) {
            n = cin.nextInt();
            BigInteger ans = BigInteger.valueOf(1);
            for (int i = 2; i <= n; i++) {
                ans = ans.multiply(BigInteger.valueOf(i));
            }
            System.out.println(ans);
        }
    }

    public static void main(String[] args) {
        Main work = new Main();
        work.solve();
    }
}
```

#### Fibonacci 数列

```java
import java.math.BigInteger;
import java.util.Scanner;

public class Main {
    void solve() {
        Scanner cin = new Scanner(System.in);
        BigInteger f1, f2, f3, f4, ans;
        while (cin.hasNext()) {
            int n = cin.nextInt();
            f1 = BigInteger.valueOf(1);
            f2 = f3 = f4 = ans = f1;
            if (n <= 4) {
                System.out.println("1");
                continue;
            }
            for (int j = 5; j <= n; j++) {
                ans = f1.add(f2.add(f3.add(f4)));
                f1 = f2;
                f2 = f3;
                f3 = f4;
                f4 = ans;
            }
            System.out.println(ans);
        }
    }

    public static void main(String[] args) {
        Main work = new Main();
        work.solve();
    }
}
```

#### 高精度小数

```java
import java.math.*;
import java.util.*;

public class Main {
    void solve() {
        Scanner cin = new Scanner(System.in);
        BigDecimal a = BigDecimal.valueOf(0);
        BigDecimal b = BigDecimal.valueOf(0);
        while (cin.hasNext()) {
            a = cin.nextBigDecimal();
            b = cin.nextBigDecimal();
            System.out.println(a.add(b).stripTrailingZeros().toPlainString());
        }
    }

    public static void main(String[] args) {
        Main work = new Main();
        work.solve();
    }
}
```

## STL

### 结构体排序

```cpp
struct node {
    int x, y;
    friend bool operator <(node a, node b) {
        if(a.x == b.x) return a.y > b.y;
        return a.x > b.x;
    }
};
```

### set

#### 迭代时删除

```cpp
for (it = a.begin(); it != a.end();) {
    if (check()) {
        set<int>::iterator tit = it;
        it++;
        a.erase(tit);
    } else it++;
}
```

## 黑科技

### FastIO

```cpp
namespace fastIO {
    #define BUF_SIZE 100000
    // fread -> read
    bool IOerror = 0;
    inline char nc() {
        static char buf[BUF_SIZE], *p1 = buf + BUF_SIZE, *pend = buf + BUF_SIZE;
        if (p1 == pend) {
            p1 = buf;
            pend = buf + fread(buf, 1, BUF_SIZE, stdin);
            if (pend == p1) {
                IOerror = 1;
                return -1;
            }
        }
        return *p1++;
    }

    inline bool blank(char ch) {
        return ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t';
    }

    inline void read(int &x) {
        char ch;
        while (blank(ch = nc()));
        if (IOerror) return;
        for (x = ch - '0'; (ch = nc()) >= '0' && ch <= '9'; x = x * 10 + ch - '0');
    }
    #undef BUF_SIZE
};

while (fastIO::read(n), !fastIO::IOerror) {}
```

### int128

```cpp
__int128 a;

if (a == 0) printf("0");
int out[50];
int len = 0;
while(a > 0) {
    out[len] = a % 10;
    a /= 10;
    len++;
}
for (int i = len - 1; i >= 0; i--) printf("%d", out[i]);
printf("\n");
```

## vimrc

```vimrc
set bs=indent,eol,start cin nu ts=4 sw=4
map <F4> :w<cr>:!gedit %<cr><cr>
map <F5> :w<Cr>:!g++ %
map <F6> <F5>;./a.out<cr>
map <F7> :!vim in<cr><cr>
map <F8> <F5>;./a.out < in<cr>
map <F9> <F5>;./a.out < in > out<cr><cr>
map <F10> :!vim out<cr><cr>
```
