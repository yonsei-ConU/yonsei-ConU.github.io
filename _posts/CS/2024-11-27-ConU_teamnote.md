---
title: "파이썬 알고리즘 정리"
excerpt: "지금까지 구현한 모든 파이썬 알고리즘 구현체 집합"
categories:
    - CS
toc: false
toc_sticky: false
date: 2024-11-27
last_modified_at: 2024-12-17
---
```py
import sys
from collections import deque
input_ = sys.stdin.readline
def minput(): return map(int, input_().split())


"""
1. 수학
1-1. isprime(n: int) -> bool
    318665857834031151167461 미만의 자연수에 대해 결정적임
    miller_rabin(a, n) 함수 동반 필요
1-2. miller_rabin(a, n)
    isprime(n)함수의 helper function
    밀러-라빈 소수판별법 구현체
1-3. pollard_rho(n: int) -> int
    폴라드-로 인수분해 알고리즘 구현체
    소인수를 하나 찾아서 돌려줌
1-4. euler_phi(n: int, l: set[int]) -> int
    phi(n) := (n 이하 자연수 중 n과 서로소인 자연수의 개수)를 반환
1-5. matrix_mult(a: matrix, b: matrix, mod: int) -> matrix
    두 행렬 a, b의 곱 행렬을 mod로 나눈 나머지 계산
1-6. matrix_pow(base: matrix, exponent: int, mod: int) -> int
    base 행렬의 exponent 제곱을 mod로 나눈 나머지 계산
1-7. sieve(n: int)
    에라토스테네스의 체 구현체
    prime_check 또는 primes를 적절히 반환
1-8. fft(arr: list, inverse: bool) -> list
    conv의 helper function
    주어진 arr의 discrete fourier transform을 반환
    나머지가 mod이면 a * 2**b + 1, 원시근이 root
           mod      a     b   root
     998244353    119    23      3
     469762049      7    26      3
    2013265921     15    27     31
    2281701377     17    27      3
    3221225473      3    30      5
1-9. conv(a: list, b: list) -> list:
    두 리스트 a, b의 이산합성곱을 반환
"""


def isprime(n):
    if n <= 71:
        if n in {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71}:
            return True
        else:
            return False
    else:
        for i in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
            if not miller_rabin(i, n):
                return False
        return True


def miller_rabin(a, n):
    d = n - 1
    r = 0
    while not d & 1:
        d >>= 1
        r += 1

    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True

    for i in range(r - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return True
    return False


def pollard_rho(n):
    from random import randint
    from math import gcd
    if isprime(n):
        return n

    for i in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]:
        if not n % i:
            return i
    g = lambda x, n, r: (x ** 2 + r) % n
    d = 1
    x = randint(2, n)
    y = x
    c = randint(1, n)

    while not d - 1:
        y = g(g(y, n, c), n, c)
        x = g(x, n, c)
        t = abs(x - y)
        d = gcd(t, n)

        if d == n:
            return pollard_rho(n)
    if isprime(d):
        return d
    return pollard_rho(d)


def euler_phi(n, l: set): # calculates phi(n), l: distinct prime divisors
    from math import gcd
    nu = n
    de = 1
    if n == 1:
        return 1
    if len(l) == 1:
        return n - 1
    else:
        for i in l:
            nu *= i - 1
            de *= i
            t = gcd(nu, de)
            nu //= t
            de //= t
        return nu // de


def matrix_mult(a, b, mod):
    r = [[0] * len(b[0]) for _ in range(len(a))]
    for p in range(len(a)):
        for q in range(len(b[0])):
            for s in range(len(a[0])):
                r[p][q] += a[p][s] * b[s][q]
            r[p][q] %= mod
    return r


def matrix_pow(base, exponent, mod):
    ret = []
    for i in range(len(base)):
        lst = [0] * len(base)
        lst[i] = 1
        ret.append(lst)
    while exponent:
        if exponent & 1:
            ret = matrix_mult(ret, base, mod)
        exponent >>= 1
        base = matrix_mult(base, base, mod)
    return ret


def sieve(n):
    prime_check = [False, False] + [True] * (n - 1)
    primes = []

    for i in range(2, n + 1):
        if prime_check[i]:
            primes.append(i)
            for j in range(i * i, n + 1, i):
                prime_check[j] = False
    return tuple(primes)


def fft(arr, inverse=False):
    mod = 998244353
    root = 3
    n = len(arr)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while not (j := j ^ bit) & bit:
            bit >>= 1
        if i < j:
            arr[i], arr[j] = arr[j], arr[i]

    length = 2
    while length <= n:
        omega = pow(root, (mod - 1) // length, mod)
        if inverse:
            omega = pow(omega, -1, mod)
        for i in range(0, n, length):
            w = 1
            for j in range(length // 2):
                u = arr[i + j]
                v = arr[i + j + length // 2] * w % mod
                arr[i + j] = (u + v) % mod
                arr[i + j + length // 2] = (u - v) % mod
                w = w * omega % mod
        length *= 2

    if inverse:
        inv_n = pow(n, -1, mod)
        arr = [(x * inv_n) % mod for x in arr]
    return arr


def conv(a, b):
    mod = 998244353
    n = len(a) + len(b) - 1
    size = 1
    while size < n:
        size <<= 1
    a.extend([0] * (size - len(a)))
    b.extend([0] * (size - len(b)))
    a = fft(a)
    b = fft(b)
    c = [(a[i] * b[i]) % mod for i in range(size)]
    c = fft(c, True)
    return c[:n]


"""
2. 기하
2-1. line_segment_intersection(a: list, b: list, c: list, d: list) -> bool
    a, b, c, d는 [x, y]꼴이어야 함
    선분 ab, 선분 cd가 교차하는지를 반환
2-2. 다각형의 넓이
2-2-1. polygon_area(x: list, y: list) -> float
    i번점은 (x[i], y[i])
2-2-2. polygon_area(points: list[list[int]]) -> float
    points[i] = (x, y)
2-3. convex_hull(points: list[list[int]]) -> list[list[int]]
    Monotone Chain 알고리즘 구현체
    결과로 나오는 컨벡스헐은 반시계방향으로 정렬되어 있음
2-4. rotating_calipers(hull)
    가먼두를 구하는 회전하는 캘리퍼스 구현체
    가장 먼 두 점 사이 거리와 실제로 가장 먼 두 점을 리턴
2-5. point_in_convex_polygon(p: [x, y], points: list[list[int]]) -> bool
    O(log N) 볼록 다각형 내부의 점 판정 구현체
    다각형의 모서리나 꼭짓점 위에 판정하려는 점이 있어도 True
2-6. angle_sort(points: list[list[int]], center: [x, y]) -> list[list[int]]
    각도 정렬 구현체
    pi/2 + epsilon이 가장 작은 각
2-7. point_in_non_convex_polygon(p: [x, y], points: list[list[int]]) -> bool
    O(N) 다각형 내부의 점 판정
    좌표범위가 10**12보다 큰 경우 수정해야 할 수도 있음
    볼록, 오목 다각형 모두에 대해 사용 가능
2-8. 픽의 정리
    넓이 = 내부격자점수 + 둘레격자점수/2 - 1
    이대로 구현하면 오차나므로 적당히 이항해야됨
"""


def line_segment_intersection(a, b, c, d):
    def ccw(a, b, c):
        cross_product = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if cross_product > 0:
            return 1
        elif cross_product < 0:
            return -1
        else:
            return 0

    ab = ccw(a, b, c) * ccw(a, b, d)
    cd = ccw(c, d, a) * ccw(c, d, b)

    if ab == 0 and cd == 0:
        a, b = sorted([a, b])
        c, d = sorted([c, d])
        return not (b < c or d < a)

    return ab <= 0 and cd <= 0


def polygon_area(x, y):
    """x and y coordinates"""
    r = 0
    for i in range(len(x)-1):
        r += x[i] * y[i+1]
    r += x[-1] * y[0]
    for i in range(len(y)):
        r -= x[i] * y[i-1]
    return r / 2


def polygon_area(points):
    """list of 2tuples"""
    area = 0
    for i in range(len(points)):
        p1 = points[i]
        p2 = points[(i + 1) % len(points)]
        area += p1[0] * p2[1] - p1[1] * p2[0]
    return abs(area) / 2


def convex_hull(points):
    if len(points) < 3: return points
    ccw = lambda p1, p2, p3: p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p2[0] * p1[1] - p3[0] * p2[1] - p1[0] * p3[1]
    points.sort()
    ret_down = [points[0], points[1]]
    for p in points:
        while len(ret_down) > 1 and ccw(ret_down[-2], ret_down[-1], p) <= 0:
            ret_down.pop()
        ret_down.append(p)

    ret_up = [points[-1], points[-2]]
    for p in points[::-1]:
        while len(ret_up) > 1 and ccw(ret_up[-2], ret_up[-1], p) <= 0:
            ret_up.pop()
        ret_up.append(p)

    return ret_up[:-1] + ret_down[:-1]


def rotating_calipers(hull):
    # 최대거리, 거리가 최대인 두 점을 리턴
    n = len(hull)
    dist2 = lambda p1, p2: (p2[1] - p1[1]) ** 2 + (p2[0] - p1[0]) ** 2
    ccw = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    if n == 1:
        return 0, hull[0], hull[0]
    elif n == 2:
        return dist2(hull[0], hull[1]), hull[0], hull[1]
    elif n == 3:
        d1 = dist2(hull[0], hull[1])
        d2 = dist2(hull[0], hull[2])
        d3 = dist2(hull[1], hull[2])
        if d1 >= d2 and d1 >= d3:
            return d1, hull[0], hull[1]
        elif d2 >= d1 and d2 >= d3:
            return d2, hull[0], hull[2]
        else:
            return d3, hull[1], hull[2]

    max_dist = 0
    actual_points = []
    j = 1
    for i in range(n):
        next_i = (i + 1) % n
        while True:
            next_j = (j + 1) % n
            cross = ccw(hull[i], hull[next_i], hull[next_j]) - ccw(hull[i], hull[next_i], hull[j])
            if cross > 0:
                j = next_j
            else:
                break

        current_dist = dist2(hull[i], hull[j])
        if current_dist > max_dist:
            max_dist = current_dist
            actual_points = [hull[i], hull[j]]
    return max_dist, *actual_points


def point_in_convex_polygon(p, polygon):
    """
    O(log N)
    assuming that polygon is ccw
    """
    ccw = lambda p1, p2, p3: p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p2[0] * p1[1] - p3[0] * p2[1] - p1[0] * p3[1]
    rlwns = polygon[0]
    lo = 0
    hi = len(polygon)
    while lo + 1 < hi:
        mid = (lo + hi) >> 1
        t = ccw(rlwns, polygon[mid], p)
        if t > 0:
            lo = mid
        elif t < 0:
            hi = mid
        else:
            if ((p[0] - rlwns[0]) ** 2 + (p[1] - rlwns[1]) ** 2 <= (polygon[mid][0] - rlwns[0]) ** 2 + (polygon[mid][1] - rlwns[1]) ** 2) and ((p[0] - polygon[mid][0]) ** 2 + (p[1] - polygon[mid][1]) ** 2 <= (polygon[mid][0] - rlwns[0]) ** 2 + (polygon[mid][1] - rlwns[1]) ** 2):
                return True
            else:
                return False
    if not lo or hi == len(polygon):
        return False
    return ccw(polygon[lo], polygon[hi], p) >= 0


def angle_sort(points, center):
    from functools import cmp_to_key
    def ccw(p1, p2, p3):
        return p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p2[0] * p1[1] - p3[0] * p2[1] - p1[0] * p3[1]

    def cmp(a, b):
        if (a < center) == (b < center):
            c = ccw(center, a, b)
            if c > 0:
                return -1
            elif c < 0:
                return 1
            else:
                return 0
        elif a < b:
            return -1
        else:
            return 1

    return sorted(points, key=cmp_to_key(cmp))


def point_in_non_convex_polygon(p, polygon):
    """진짜 역대급개쓰레기"""
    sp0 = [p[0] + 0.000000000001, p[1] + 0.000000000001]
    sp1 = [p[0] + 1, p[1] + 10 ** 12]
    tmp = 0
    for i in range(len(polygon) - 1):
        if p == polygon[i]:
            return True
        if line_segment_intersection(sp0, sp1, polygon[i], polygon[i + 1]):
            tmp += 1
    if p == polygon[-1]:
        return True
    if line_segment_intersection(sp0, sp1, polygon[-1], polygon[0]):
        tmp += 1
    return tmp & 1 == 1


"""
3. 자료 구조
3-1. UnionFind
    분리 집합 구현체
    UnionFind.__init__(x): 크기 x의 분리 집합 초기화
    UnionFind.find(x): 경로 압축 최적화 적용, x정점의 루트노드 반환
    UnionFind.union(x, y): x, y정점을 합침
3-2. segtree
    비재귀 세그먼트 트리 구현체
    연산에 결합법칙이 성립하지 않아도 작동함
    segtree.__init__(arr, func, identity): 초기 배열이 arr, 합치는 함수가 func, func의 항등원이 identity인 세그먼트 트리를 만듦
    segtree.update(idx, val): idx번째 인덱스 값을 val로 바꿈
    segtree.query(l, r): l번째부터 r번째까지 인덱스에 func를 적용한 결과를 반환
3-3. lazy_segtree
    Lazy Propagation 구현체
    아직 덧셈만 구현함
    사용법은 segtree와 동일
"""


class UnionFind:
    def __init__(self, size):
        self.parent = [i for i in range(size)]
        self.rank = [0] * size

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)
        if root_x != root_y:
            if self.rank[root_x] > self.rank[root_y]:
                self.parent[root_y] = root_x
            elif self.rank[root_x] < self.rank[root_y]:
                self.parent[root_x] = root_y
            else:
                self.parent[root_y] = root_x
                self.rank[root_x] += 1


class segtree:
    """ConU's non-recursive uniform segment tree implementation"""
    def __init__(self, arr, func, identity):
        i = 1
        while i < len(arr): i <<= 1
        self.n = i
        self.tree = [identity for _ in range(self.n << 1)]
        self.func = func
        self.identity = identity
        for i in range(len(arr)):
            self.tree[self.n + i] = arr[i]
        for i in range(self.n - 1, 0, -1):
            self.tree[i] = func(self.tree[i << 1], self.tree[(i << 1) | 1])

    def update(self, idx, val):
        idx += self.n
        self.tree[idx] = val
        while idx > 1:
            idx >>= 1
            self.tree[idx] = self.func(self.tree[idx << 1], self.tree[(idx << 1) | 1])

    def query(self, l, r):
        ret_left = self.identity
        ret_right = self.identity
        l += self.n
        r += self.n
        while l <= r:
            if l & 1:
                ret_left = self.func(ret_left, self.tree[l])
                l += 1
            if not r & 1:
                ret_right = self.func(self.tree[r], ret_right)
                r -= 1
            l >>= 1
            r >>= 1
        return self.func(ret_left, ret_right)


class lazy_segtree:
    """only lazy sum seg"""
    def __init__(self, arr):
        self.n = 1
        while self.n < len(arr):
            self.n <<= 1
        self.tree = [0] * (2 * self.n)
        self.lazy = [0] * (2 * self.n)
        for i in range(len(arr)):
            self.tree[self.n + i] = arr[i]
        for i in range(self.n - 1, 0, -1):
            self.tree[i] = self.tree[2 * i] + self.tree[2 * i + 1]

    def propagate(self, node, left, right):
        if self.lazy[node] != 0:
            self.tree[node] += (right - left + 1) * self.lazy[node]
            if left != right:
                for child in [2 * node, 2 * node + 1]:
                    self.lazy[child] += self.lazy[node]
            self.lazy[node] = 0

    def update(self, l, r, add=0, node=1, left=0, right=None):
        if right is None:
            right = self.n - 1
        self.propagate(node, left, right)
        if r < left or right < l:
            return
        if l <= left and right <= r:
            self.lazy[node] += add
            self.propagate(node, left, right)
            return
        mid = (left + right) // 2
        self.update(l, r, add, 2 * node, left, mid)
        self.update(l, r, add, 2 * node + 1, mid + 1, right)
        self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]

    def query(self, l, r, node=1, left=0, right=None):
        if right is None:
            right = self.n - 1
        self.propagate(node, left, right)
        if r < left or right < l:
            return 0
        if l <= left and right <= r:
            return self.tree[node]
        mid = (left + right) // 2
        p1 = self.query(l, r, 2 * node, left, mid)
        p2 = self.query(l, r, 2 * node + 1, mid + 1, right)
        return p1 + p2


"""
4. DP
4-1. LIS_len(arr)
    가장 긴 증가하는 부분 수열 길이를 리턴 O(NlogN)
4-2. LIS(arr)
    가장 긴 증가하는 부분 수열을 아무거나 하나 리턴
4-3. LCS_len(str1, str2)
    두 문자열의 longest increasing subsequence의 길이를 구함
4-4. LCS(str1, str2)
    두 문자열의 LCS를 아무거나 리턴
    세 문자열의 LCS를 원할 때는 LCS(LCS(str, str2), str3)으로 구하면 안 됨
"""


def LIS_len(arr):
    temp = [-float('inf')]
    for element in arr:
        lo = -1
        hi = len(temp)
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if temp[mid] < element:
                lo = mid
            else:
                hi = mid
        if hi == len(temp):
            temp.append(element)
        else:
            temp[hi] = element
    return len(temp) - 1


def LIS(arr):
    temp = [-float('inf')]
    idx = []
    for element in arr:
        lo = -1
        hi = len(temp)
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if temp[mid] < element:
                lo = mid
            else:
                hi = mid
        if hi == len(temp):
            temp.append(element)
        else:
            temp[hi] = element
        idx.append(hi)
    result = []
    for i in range(len(arr)-1, -1, -1):
        if idx[i] == len(temp) - 1 - len(result):
            result.append(arr[i])
    return result[::-1]


def LCS_len(str1, str2):
    len1 = len(str1)
    len2 = len(str2)
    dp = [[0 for i in range(len2+1)] for j in range(len1+1)]

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    return dp


def LCS(str1, str2):
    len1 = len(str1)
    len2 = len(str2)
    dp = [[0 for i in range(len2+1)] for j in range(len1+1)]

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    res = []
    i = len1
    j = len2
    while i > 0 and j > 0:
        t = dp[i][j]
        if dp[i - 1][j] == t:
            i -= 1
        elif dp[i][j - 1] == t:
            j -= 1
        else:
            res.append(str1[i-1])
            i -= 1; j -= 1
    return ''.join(res)[::-1]


"""
5. 그래프
5-1. 트리 관련 알고리즘
5-1-1. kruskal(v, edges)
    크루스칼 알고리즘 구현체
    v개의 정점을 가진 그래프, edges에는 (가중치, 시작, 끝)이 들어감
    유니온 파인드 자료구조 필요
    알고리즘 종료 이후 mst_weight에는 최소 스패닝 트리의 가중치가, mst에는 아무 mst가 들어가 있음
5-1-2. LCA_preprocess(g, root)
    sparse table을 이용한 LCA 구현 전처리 함수
    D는 루트와의 거리, d는 깊이, table은 LCA를 구하기 위한 sparse table에 해당됨
5-1-3. LCA_query(u, v, d, table)
    LCA_preprocess 함수가 미리 한 번 작동했어야 함
    정점 u, v의 최소 공통 조상 번호를 리턴
5-1-4. ETT(g, root)
    루트 정점 번호가 root일 때 오일러 경로 테크닉 구현체
    i번 정점을 루트로 하는 서브트리는 disc[i]번부터 esc[i]번까지
5-1-5. heavy_light_decomposition(graph, root)
    heavy-light 분할 구현체
    i번 정점을 루트로 하는 서브트리는 disc[i]번부터 esc[i]번까지
    나를 루트로 하는 서브트리의 크기는 size
    자기 체인의 가장 위 정점은 top
5-2. 최단 경로 알고리즘
5-2-1. dijkstra(g, st)
    다익스트라 알고리즘 구현체. g[cur] = (nxt, dist)
    시작정점 번호가 st일 때 나머지 모든 정점까지의 최단경로를 리턴
5-2-2. dijkstra_with_path(g, st)
    역추적이 포함된 다익스트라
    last[i] = st에서 i로 가는 최단경로에서 i번 직전에 거쳐가는 정점 번호
5-2-3. floyd_warshall(g)
    플로이드-워셜 알고리즘 구현체. g[cur][nxt] = dist
    모든 정점에서 모든 정점까지 최단경로의 길이
5-2-4. bellman_ford(g, st)
    벨만-포드 알고리즘 구현체. g[cur] = (nxt, dist)
    음수 사이클이 존재할 경우 False를 리턴함
5-2-5. SPFA(g, st)
    SPFA 구현체
    사용법은 벨만-포드와 동일
5-3. 네트워크 플로우
5-3-1. FlowEdge
    플로우 그래프 간선 클래스로, 다른 플로우 관련 함수들에는 이 간선 클래스가 사용되어야 함
5-3-2. add_edge(g, start, end, capacity)
    플로우 그래프 g에 start에서 end로 가는 용량이 capacity인 간선을 추가
5-3-3. dinitz(g, source, sink)
    최대 유량 Dinitz 알고리즘 구현체
5-4. 기타 그래프 알고리즘
5-4-1. topological_sorting(graph, indegree)
    위상 정렬 구현체
    그래프와, 각 정점의 indegree를 저장한 리스트를 받아서 위상정렬된 정점 순서를 리턴
5-4-2. DFS
5-4-2-1. dfs_recursive(begin, connect)
    재귀 DFS 구현체, A~E 함수 위치에 주의
5-4-2-2. dfs_nonrecursive(begin, connect)
    비재귀 DFS 구현체, 재귀 DFS에서 A~E위치에 들어갈 함수를 넣으면 됨
5-4-3. tarjan(g)
    SCC 타잔 알고리즘 구현체
    강한 연결 요소가 위상정렬된 상태로 리턴됨
5-4-4. two_sat(N, clauses, trace=False)
    변수 N개짜리 2-sat의 만족 가능 여부를 리턴
    trace값이 참이면, 모순이 일어나지 않을 때 가능한 실제 해를 리턴
    tarjan 함수 필요
    clauses에는 11280번 문제에서 주어지는 것과 같은 형식으로 입력
5-4-5. find_articulation_point(g)
    단절점 알고리즘 구현체
    모든 단절점을 리턴
5-4-6. find_bridge(g)
    단절선 알고리즘 구현체
    모든 단절선을 리턴
5-4-7. hopcroft_karp(adj, n, m)
    호프크로프트-카프 이분 매칭 알고리즘 구현체
    adj는 그래프, n은 출발 정점 개수, m은 도착 정점 개수
"""


def kruskal(v, edges):
    edges.sort()
    uf = UnionFind(v + 1)
    mst = set()
    mst_weight = 0
    for edge in edges:
        weight, s, e = edge
        if uf.find(s) == uf.find(e):
            continue
        else:
            mst.add(edge)
            mst_weight += weight
            uf.union(s, e)
    return mst_weight


def LCA_preprocess(g, root):
    v = len(g)
    D = [0] * v
    d = [0] * v
    table = [[0] * v for _ in range(v.bit_length())]

    def LCA_dfs(cur, parent):
        table[0][cur] = parent
        for nxt, weight in g[cur]:
            if nxt == parent:
                continue
            D[nxt] = D[cur] + weight
            d[nxt] = d[cur] + 1
            LCA_dfs(nxt, cur)

    LCA_dfs(root, -1)

    for i in range(1, v.bit_length()):
        for j in range(v):
            table[i][j] = table[i-1][table[i-1][j]]

    return D, d, table  # distance, depth, table


def LCA_query(u, v, d, table):
    if d[v] > d[u]:
        u, v = v, u
    x = d[u] - d[v]
    for i in range(x.bit_length()):
        if x & 1:
            u = table[i][u]
        x >>= 1
    if u == v:
        return u
    for j in range(len(table) - 1, -1, -1):
        if table[j][u] != table[j][v]:
            u = table[j][u]
            v = table[j][v]
    return table[0][v]


def ETT(g, root):
    time = -1
    disc = [-1] * len(g)
    esc = [-1] * len(g)

    def ETT_process(cur, parent):
        nonlocal time
        time += 1
        disc[cur] = time
        for nxt in g[cur]:
            if nxt == parent:
                continue
            ETT_process(nxt, cur)
        esc[cur] = time
    ETT_process(root, -1)
    return disc, esc


def heavy_light_decomposition(graph, root):
    v = len(graph)
    sizes = [0] * v
    depth = [0] * v
    parent = [root] * v
    disc = [-1] * v
    esc = [-1] * v
    top = [root] * v
    time = -1
    g = [[] for _ in range(v)]

    def get_g(cur, par):
        for nxt in graph[cur]:
            if nxt == par:
                continue
            g[cur].append(nxt)
            get_g(nxt, cur)

    def decompose(cur):
        sizes[cur] = 1
        for nxt in g[cur]:
            depth[nxt] = depth[cur] + 1
            parent[nxt] = cur
            decompose(nxt)
            sizes[cur] += sizes[nxt]
        g[cur].sort(key=lambda x: sizes[x], reverse=True)

    def hld_ett(cur):
        nonlocal time
        time += 1
        disc[cur] = time
        for nxt in g[cur]:
            top[nxt] = top[cur] if nxt == g[cur][0] else nxt
            hld_ett(nxt)
        esc[cur] = time

    get_g(root, -1)
    decompose(root)
    hld_ett(root)

    return sizes, depth, parent, disc, esc, top


def dijkstra(g, st):
    '''
    g: graph, g[start] = (end, dist)
    st: start node
    '''
    import heapq
    distances = [float('inf')] * len(g)
    distances[st] = 0
    heap = []
    heapq.heappush(heap, (0, st))

    while heap:
        dist, cur = heapq.heappop(heap)
        if distances[cur] < dist:
            continue
        for nextnum, nextdist in g[cur]:
            t = dist + nextdist
            if distances[nextnum] > t:
                distances[nextnum] = t
                heapq.heappush(heap, (t, nextnum))

    return distances


def dijkstra_with_path(g, st):
    import heapq
    v = len(g)
    distances = [float('inf')] * v
    last = [0] * v
    distances[st] = 0
    heap = []
    heapq.heappush(heap, (0, st))

    while heap:
        dist, cur = heapq.heappop(heap)
        if distances[cur] < dist:
            continue
        for nextnum, nextdist in g[cur]:
            t = dist + nextdist
            if distances[nextnum] > t:
                distances[nextnum] = t
                last[nextnum] = cur
                heapq.heappush(heap, (t, nextnum))

    return distances, last


def floyd_warshall(g):
    """
    g: graph, g[start][end] = dist
    return: distance matrix
    """
    n = len(g)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if i != j:
                    g[i][j] = min(g[i][k] + g[k][j], g[i][j])
    return g


def bellman_ford(g, st):
    """
    g: graph, g[cur] = (nxt, dist)
    st: the vertex to start
    return: if there exists a negative cycle, False
            if there does not, distances list
    """
    v = len(g)
    distances = [float('inf')] * v
    distances[st] = 0

    for i in range(v - 1):
        for cur in range(v):
            for nxt, dist in g[cur]:
                if distances[cur] != float('inf') and distances[cur] + dist < distances[nxt]:
                    distances[nxt] = distances[cur] + dist

    for cur in range(v):
        for nxt, dist in g[cur]:
            if distances[cur] != float('inf') and distances[cur] + dist < distances[nxt]:
                return False

    return distances


def SPFA(g, st):
    """
    g: graph, g[cur] = (nxt, dist)
    st: the vertex to start
    return: if there exists a negative cycle, False
            if there does not, distances list
    """
    from collections import deque
    v = len(g)
    distances = [float('inf')] * v
    distances[st] = 0

    q = deque([st])
    in_queue = [False] * v
    in_queue[st] = True
    cnt = [0] * v
    cnt[st] = 1

    while q:
        cur = q.popleft()
        in_queue[cur] = False
        for nxt, dist in g[cur]:
            if distances[cur] + dist < distances[nxt]:
                distances[nxt] = distances[cur] + dist
                if not in_queue[nxt]:
                    in_queue[nxt] = True
                    q.append(nxt)
                    cnt[nxt] += 1
                    if cnt[nxt] >= v:
                        return False

    return distances


class FlowEdge:
    def __init__(self, start, end, capacity):
        self.start = start
        self.end = end
        self.capacity = capacity
        self.flow = 0
        self.reverse = None

    def __repr__(self):
        return f"({self.end}, capacity: {self.capacity}, flow: {self.flow})"


def add_edge(g, start, end, capacity):
    forward = FlowEdge(start, end, capacity)
    backward = FlowEdge(end, start, 0)
    forward.reverse = backward
    backward.reverse = forward
    g[start].append(forward)
    g[end].append(backward)


def dinitz(g, source, sink):
    N = len(g)
    level = [-1] * N
    ret = 0

    def dinitz_bfs():
        nonlocal level
        level = [-1] * N
        q = deque([source])
        level[source] = 0
        while q:
            cur = q.popleft()
            for edge in g[cur]:
                if level[edge.end] == -1 and edge.flow < edge.capacity:
                    level[edge.end] = level[cur] + 1
                    q.append(edge.end)
        return level[sink] != -1

    def dinitz_dfs(cur, flow):
        if cur == sink:
            return flow
        while ptr[cur] < len(g[cur]):
            edge = g[cur][ptr[cur]]
            if level[edge.end] == level[cur] + 1 and edge.flow < edge.capacity:
                pushed = dinitz_dfs(edge.end, min(flow, edge.capacity - edge.flow))
                if pushed > 0:
                    edge.flow += pushed
                    edge.reverse.flow -= pushed
                    return pushed
            ptr[cur] += 1
        return 0

    while dinitz_bfs():
        ptr = [0] * N
        while True:
            pushed = dinitz_dfs(source, float('inf'))
            if pushed == 0:
                break
            ret += pushed

    return ret


def topological_sort(graph, indegree):
    zero_indegree = deque([i for i in range(len(graph)) if not indegree[i]])
    free = [True] * len(indegree)

    result = []
    while zero_indegree:
        node = zero_indegree.popleft()
        result.append(node)
        free[node] = False
        for neighbor in graph[node]:
            indegree[neighbor] -= 1
            if indegree[neighbor] == 0:
                zero_indegree.append(neighbor)

    assert len(result) == len(graph)

    return result + [i for i in range(len(free)) if i and free[i]]


def dfs_recursive(u, connect, visited):
    """
    u: start vertex
    connect: graph, adjacency list
    visited: visited
    A, B, C, D, E: any function you want
    """
    A = B = C = D = E = lambda: 1
    A(u)
    visited[u] = True
    for v in connect[u]:
        if not visited[v]:
            B(u, v)
            dfs_recursive(v)
            C(u, v)
        else:
            E(u, v)
    D(u)


def dfs_nonrecursive(begin, connect):
    """
    begin: start vertex
    connect: adjacency list
    A, B, C, D, E: any function you want, corresponds to `dfs_recursive`
    """
    A = B = C = D = E = lambda: 1
    visited = [False] * len(connect)
    stack = [(begin, -1)]  # 정점 번호, 인접 리스트의 이웃 정점 인덱스
    while stack:
        u, i = stack.pop()
        visited[u] = True
        if i < 0:
            A(u)
        else:
            v = connect[u][i]
            C(u, v)
        i += 1
        while i < len(connect[u]):
            v = connect[u][i]
            if visited[v]:
                E(u, v)
                i += 1
                continue
            B(u, v)
            stack.append((u, i))
            stack.append((v, -1))
            break
        else:
            D(u)


def tarjan(g):
    v = len(g)
    disc = [-1] * v
    low = [-1] * v
    stack = []
    in_stack = [False] * v
    t = 0
    ret = []

    def dfs(cur):
        nonlocal t
        disc[cur] = t
        low[cur] = t
        t += 1
        stack.append(cur)
        in_stack[cur] = True

        for nxt in g[cur]:
            if disc[nxt] == -1:
                dfs(nxt)
                low[cur] = min(low[cur], low[nxt])
            elif in_stack[nxt]:
                low[cur] = min(low[cur], disc[nxt])

        if disc[cur] == low[cur]:
            scc = []
            while 1:
                node = stack.pop()
                scc.append(node)
                in_stack[node] = False
                if cur == node:
                    break
            ret.append(scc)

    for i in range(v):
        if disc[i] == -1:
            dfs(i)

    return ret


def two_sat(N, clauses, trace=False):
    """
    변수개수 N
    clauses는 list[pii], 각 pii에 들어있는 정수가 i일때, i > 0이면 i번변수가 참, i < 0이면 i번변수가 거짓
    모순이 일어난다면 False
    그렇지 않다면
        trace가 False라면 1
        trace가 True라면 각 명제별로 참이면 1, 거짓이면 0이 되는 해 리스트
    """
    g = [[] for _ in range(N << 1)]
    for a, b in clauses:
        if a < 0:
            a = N - a
        if b < 0:
            b = N - b
        a -= 1; b -= 1
        g[(a + N) % (2 * N)].append(b)
        g[(b + N) % (2 * N)].append(a)
    scc = tarjan(g)
    scc_rev = [0] * (N << 1)
    for i in range(len(scc)):
        for c in scc[i]:
            scc_rev[c] = i
    for i in range(N):
        if scc_rev[i] == scc_rev[i + N]:
            return 0
    if not trace:
        return 1
    ret = [-1] * N
    for component in scc:
        for c in component:
            idx = c if c < N else c - N
            if ret[idx] == -1:
                ret[idx] = c < N
    return ret


def find_articulation_point(g):
    time = -1
    disc = [-1] * len(g)
    low = [-1] * len(g)
    ret = [False] * len(g)

    def dfs(cur, parent):
        nonlocal time
        time += 1
        disc[cur] = time
        low[cur] = time
        children_count = 0
        for nxt in g[cur]:
            if disc[nxt] == -1:
                children_count += 1
                dfs(nxt, cur)
                low[cur] = min(low[cur], low[nxt])
                if parent != -1 and low[nxt] >= disc[cur]:
                    ret[cur] = True
            elif nxt != parent:
                low[cur] = min(low[cur], disc[nxt])
        if parent == -1 and children_count > 1:
            ret[cur] = True

    for root in range(len(g)):
        if disc[root] == -1:
            dfs(root, -1)

    return ret


def find_bridge(g):
    time = -1
    disc = [-1] * len(g)
    low = [-1] * len(g)
    ret = []

    def dfs(cur, parent):
        nonlocal time
        time += 1
        disc[cur] = time
        low[cur] = time
        for nxt in g[cur]:
            if disc[nxt] == -1:
                dfs(nxt, cur)
                low[cur] = min(low[cur], low[nxt])
                if low[nxt] > disc[cur]:
                    c, n = cur, nxt
                    if c > n: c, n = n, c
                    ret.append((c, n))
            elif nxt != parent:
                low[cur] = min(low[cur], disc[nxt])

    for root in range(len(g)):
        if disc[root] == -1:
            dfs(root, -1)

    return ret


def hopcroft_karp(adj, n, m):
    from collections import deque

    g_match = [-1] * n
    h_match = [-1] * m
    dist = [0] * n
    inf = float('inf')

    def hk_bfs():
        q = deque()
        for u in range(n):
            if g_match[u] == -1:
                dist[u] = 0
                q.append(u)
            else:
                dist[u] = inf

        ret = inf
        while q:
            u = q.popleft()
            if dist[u] >= ret:
                continue
            for v in adj[u]:
                if h_match[v] == -1:
                    ret = dist[u] + 1
                elif dist[h_match[v]] == inf:
                    dist[h_match[v]] = dist[u] + 1
                    q.append(h_match[v])

        return ret != inf

    def hk_dfs(u):
        for v in adj[u]:
            if h_match[v] == -1 or (dist[h_match[v]] == dist[u] + 1 and hk_dfs(h_match[v])):
                g_match[u] = v
                h_match[v] = u
                return True
        dist[u] = inf
        return False

    ret = 0
    while hk_bfs():
        for u in range(n):
            if g_match[u] == -1:
                ret += hk_dfs(u)

    return ret


"""
6. 문자열
6-1. knuth_morris_pratt(s1, s2)
    KMP 알고리즘 구현체
    s1 문자열에서 s2 문자열을 검색
    s1 문자열에서 s2 문자열의 시작점이 나타나는 인덱스의 리스트를 리턴
    중간에 fail함수를 리턴하면 실패함수만 얻을 수 있음
6-2. TrieNode
    트라이 자료구조를 위한 노드 구현체
    이것 자체로는 쓸 일이 없고 trie.insert할 때 알아서 사용됨
6-3. Trie
    트라이 자료구조 구현체
    trie.insert(word): word 단어를 트라이에 추가. TrieNode 클래스 사용
    trie.search(word): word 단어를 트라이에서 검색. 이 부분 알잘딱하게 바꾸면 이것저것 할 수 있음
6-4. manacher(s)
    매내처 알고리즘 구현체
    중간에 더미 문자열을 삽입하고, 각 문자별로 그 문자열을 중심으로 하는 가장 긴 팰린드롬 부분 문자열의 길이를 리턴
6-5. z(s)
    Z 알고리즘 구현체
    z[i]: i번째 인덱스에서 시작하는 접미사와 전체 문자의 가장 긴 공통 접두사의 길이
6-6. suffix_array(S)
    Manber-Myers 알고리즘 구현체
    로그제곱 시간에 접미사 배열을 구함
    radix sort로 O(NlogN)이 될 수 있긴 한데 언제 바꾸냐
6-7. kasai(S, sa)
    문자열 S와 미리 계산된 접미사 배열 sa가 주어질 때 LCP배열을 구함
    보통 세그먼트 트리 RmQ로 구해야 함
6-8. AhoCorasickTrieNode
    아호-코라식 알고리즘을 위한 트라이 노드 클래스
6-9. AhoCorasickTrie
    아호-코라식 알고리즘을 위한 트라이 구현체
    collections.deque 임포트 필요
    AhoCorasickTrie.insert(word): word를 문자열 집합에 추가
    AhoCorasickTrie.build_aho_corasick(): insert를 다 한 뒤 BFS로 아호-코라식 알고리즘 전처리
    AhoCorasickTrie.search(text): text에서 문자열 집합에 포함된 인덱스들을 리턴
"""


def knuth_morris_pratt(s1, s2):
    fail = [0] * len(s2)
    j = 0

    for i in range(1, len(s2)):

        while j > 0 and s2[i] != s2[j]:
            j = fail[j-1]

        if s2[i] == s2[j]:
            j += 1
            fail[i] = j

    result = []
    j = 0

    for i in range(len(s1)):

        while j > 0 and s1[i] != s2[j]:
            j = fail[j-1]

        if s1[i] == s2[j]:
            if j + 1 == len(s2):
                result.append(i + 2 - len(s2))
                j = fail[j]
            else:
                j += 1

    return result


class TrieNode:
    def __init__(self, s=None):
        self.children = {}
        self.s = s

    def __str__(self):
        return str(self.s)


class Trie:
    """Do not forget TrieNode"""
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word):
        cur = self.root
        for c in word:
            if c not in cur.children:
                cur.children[c] = TrieNode(c)
            cur = cur.children[c]

    def search(self, word):
        cur = self.root
        for c in word:
            if c not in cur.children:
                return False
            cur = cur.children[c]
        return True


def manacher(s):
    dummy = '#'
    S = []
    for i in range(len(s)): S += [dummy, s[i]]
    S.append(dummy)

    ret = [0] * len(S)
    j = r = 0
    for i in range(1, len(S) - 1):
        mirror = (j << 1) - i
        if i < r: ret[i] = min(ret[mirror], r - i)
        while i + ret[i] + 1 < len(S) and i - ret[i] - 1 >= 0 and S[i + ret[i] + 1] == S[i - ret[i] - 1]: ret[i] += 1
        if i + ret[i] > r: j = i; r = i + ret[i]
    return ret


def z(s):
    ret = [0] * len(s)
    l = r = 0
    ret[0] = len(s)
    for i in range(1, len(s)):
        if i > r:
            l = r = i
            while r < len(s) and s[r - l] == s[r]: r += 1
            r -= 1
            ret[i] = r - l + 1
        else:
            if ret[i - l] < r - i + 1: ret[i] = ret[i - l]
            else:
                l = i
                while r < len(s) and s[r - l] == s[r]: r += 1
                r -= 1
                ret[i] = r - l + 1
    return ret


def suffix_array(S):
    ret = list(range(len(S)))
    rank = [ord(i) for i in S]
    tmp = [0] * len(S)
    k = 1
    while k < len(S):
        ret.sort(key=lambda i: (rank[i], rank[i + k] if i + k < len(S) else -1))
        tmp[ret[0]] = 0
        for i in range(1, len(S)):
            prev = ret[i - 1]
            cur = ret[i]
            tmp[cur] = tmp[prev] + ((rank[cur], rank[cur + k] if cur + k < len(S) else -1) != (rank[prev], rank[prev + k] if prev + k < len(S) else -1))
        tmp, rank = rank, tmp
        k <<= 1
    return ret


def kasai(S, sa):
    sa_rev = [0] * len(S)
    for i in range(len(S)): sa_rev[sa[i]] = i
    ret = [0] * len(S)
    k = 0
    for i in range(len(S)):
        if sa_rev[i] == len(S) - 1:
            k = 0
            continue
        j = sa[sa_rev[i] + 1]
        while i + k < len(S) and j + k < len(S) and S[i + k] == S[j + k]: k += 1
        ret[sa_rev[i] + 1] = k
        k -= bool(k)
    return ret


class AhoCorasickTrieNode:
    def __init__(self, s=None):
        self.children = {}
        self.s = s
        self.fail = None
        self.output = []

    def __str__(self):
        return str(self.s)


class AhoCorasickTrie:
    """do not forget to import deque"""
    def __init__(self):
        self.root = AhoCorasickTrieNode()

    def insert(self, word):
        cur = self.root
        for c in word:
            if c not in cur.children:
                cur.children[c] = AhoCorasickTrieNode(c)
            cur = cur.children[c]
        cur.output.append(word)

    def build_aho_corasick(self):
        q = deque()
        for child in self.root.children.values():
            q.append(child)
            child.fail = self.root
        while q:
            cur = q.popleft()
            for c, nxt in cur.children.items():
                fail_node = cur.fail
                while fail_node and c not in fail_node.children:
                    fail_node = fail_node.fail
                nxt.fail = fail_node.children[c] if fail_node else self.root
                if nxt.fail:
                    nxt.output += nxt.fail.output
                q.append(nxt)

    def search(self, text):
        cur = self.root
        matches = []
        for i, c in enumerate(text):
            while c not in cur.children and cur != self.root:
                cur = cur.fail
            if c in cur.children:
                cur = cur.children[c]
            else:
                cur = self.root
            if cur.output:
                for pattern in cur.output:
                    matches.append((i - len(pattern) + 1, pattern))
        return matches


"""
7. 기타 알고리즘
7-1. find_cycle(f, x0)
    플로이드의 토끼와 거북이 알고리즘 구현체
    a_(n+1) = f(a_n)이며 a_1 = x0일 때 사이클을 찾아줌
    사이클 입장까지의 길이가 mu, 한 사이클의 길이가 lam
7-2. permutation_cycle_decomposition(l, permutation)
    순열 사이클 분할 알고리즘 구현체
    길이 l의 순열 permutation이 주어질 때 순열 사이클 분할을 진행
7-3. coordinate_compression(lst)
    좌표 압축 구현체
    rank 또는 compressed 모두 리턴 가능
7-4. bootstrap
    재귀함수 펴는 무언가
    재귀함수에서, @bootstrap을 붙인 다음 재귀 호출 시 yield dfs()처럼 사용
    또한, return ret 대신 yield ret 사용, 맨 마지막에 yield 써주기
"""


def find_cycle(f, x0):
    tortoise = f(x0)
    hare = f(f(x0))
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(f(hare))

    mu = 0
    tortoise = x0
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(hare)
        mu += 1

    lam = 1
    hare = f(tortoise)
    while tortoise != hare:
        hare = f(hare)
        lam += 1

    return lam, mu


def permutation_cycle_decomposition(l, permutation):
    """1-based!!!"""
    processed = [False] * l
    cycles = []
    for i in range(l):
        if processed[i]:
            continue
        cycle = []
        pointer = i
        while not processed[pointer]:
            processed[pointer] = True
            cycle.append(pointer + 1)
            pointer = permutation[pointer] - 1
        if cycle:
            cycles.append(cycle)

    return cycles


def coordinate_compression(lst):
    """
    DO NOT USE THIS IN CODEFORCES
    return anything you want
    """
    distinct = sorted(set(lst))
    rank = {distinct[i]: i for i in range(len(distinct))}
    compressed = [rank[e] for e in lst]
    return compressed


def bootstrap(f, stack=[]):
    from types import GeneratorType
    def wrappedfunc(*args, **kwargs):
        if stack:
            return f(*args, **kwargs)
        else:
            to = f(*args, **kwargs)
            while True:
                if type(to) is GeneratorType:
                    stack.append(to)
                    to = next(to)
                else:
                    stack.pop()
                    if not stack:
                        break
                    to = stack[-1].send(to)
            return to

    return wrappedfunc

```
