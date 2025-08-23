---
title: "스트레스 테스팅의 개념"
excerpt: "알고리즘 뿐만아니라 다양한 분야에서 쓰이는 디버깅 방법에 대해 알아봅시다!"
categories:
    - CS
toc: true
toc_sticky: false
date: 2025-01-30
last_modified_at: 2025-01-30
---

# 개요
문제를 풀다가 물론 한 번에 맞는다면 기분이 좋지만 흔히 말하는 '맞왜틀'을 당할 때도 가끔 있습니다.  
논리상의 오류를 발견하고 고칠 수 있는 경우라면 다행이지만, 가끔은 어디에서 틀렸는지 도저히 모르겠을 때가 있죠.  
그럴 때 반례라도 알면 참 좋을 것 같아 아무 테스트 케이스나 넣고 돌려 보면 어떨까요?  
이 아이디어를 구체화하면 상당히 강력한 디버깅 방법인 '스트레스 테스팅'을 할 수 있게 됩니다.  

# 스트레스 테스팅의 개념
스트레스 테스팅이란 간단히 "느린 코드와 자기의 풀이를 비교해서 랜덤한 반례를 만들어 내는" 방법이라고 할 수 있어요.  
그렇기 때문에 스트레스 테스팅에는 크게 두 가지의 힘수를 넣고 돌려요.  
그 중 하나는 solve 함수로, 지금 반례를 찾고 싶은 코드를 집어 넣으면 돼요.  
나머지 함수는 naive 함수로, 느리지만 정답을 확정적으로 내 놓는 코드를 짜면 돼요.  
예를 들어 O(N!) 시간 복잡도로 브루트포스 하는 것처럼, 시간 복잡도와 상관없이 무식한 방법으로 답을 어떻게든 찾는 코드를 짜시면 되는 거죠.  

# 사용 예제
```py
import sys
from random import randint
from algorithms import segtree, sieve
input_ = sys.stdin.readline
def minput(): return map(int, input_().split())
def mult(a, b): return a * b % MOD


primes = sieve(29)
MOD = 998244353


def right(data):
    st = segtree(primes, mult, 1)
    st_original = segtree(primes, mult, 1)
    ret = []

    for Q, A, B in data:
        if not Q:
            curA = st.query(A, A)
            curB = st.query(B, B)
            st.update(A, curB)
            st.update(B, curA)
        else:
            t1 = st.query(A, B)
            t2 = st_original.query(A, B)
            ret.append('YNEOS'[t1 != t2::2])

    return ret


def wrong(data, N):
    st = segtree(list(range(1, N + 1)), mult, 1)
    st_original = segtree(list(range(1, N + 1)), mult, 1)
    ret = []

    for Q, A, B in data:
        if not Q:
            curA = st.query(A, A)
            curB = st.query(B, B)
            st.update(A, curB)
            st.update(B, curA)
        else:
            t1 = st.query(A, B)
            t2 = st_original.query(A, B)
            ret.append('YNEOS'[t1 != t2::2])

    return ret


def generate():
    N = randint(3, 6)
    K = randint(3, 6)
    data = []
    for i in range(K):
        a = randint(0, 1)
        b = randint(0, N - 2)
        c = randint(b + 1, N - 1)
        data.append([a, b, c])
    return N, data


for i in range(1, 10001):
    N, data = generate()
    r = right(data)
    w = wrong(data, N)
    bsn = '\n'
    if r != w:
        print(f"Test case #{i} failed")
        print(f"Input {N} {data}")
        print(f"Expected: {r}")
        print(f"Received: {w}")
        break
else:
    print(f"10000 Test case passed")

```
위 코드는 제가 스트레스 테스팅할 때 쓰는 템플릿이에요.  
하나하나 살펴봅시다!  
일단 가져온건 디버깅하려고 짰던 코드는 아니고 제가 짰던 다른 코드를 저격하기 위해서 만들었던 스트레스이긴 합니다 ㅜㅜ  
어쨌든 right은 제코드고 wrong도 제코드긴 한데 저격하려는 코드입니다!  
여기서는 그냥 right이 브루트포스 코드고 wrong이 제가 짠 풀이라고 합시다.  
뭐 나중에 스트레스 돌릴 적절한 예시가 있다면 그때 바꿔 볼게요  
적절한 형태의 랜덤 generate() 함수를 만들어 준 다음에, 원하는 만큼 반복해서 실행해주면 됩니다!  
  
그리고 나면 아래랑 비슷한 출력을 얻을 수 있게 돼요.
```
Test case #1506 failed
Input 6 [[0, 0, 3], [0, 1, 5], [0, 0, 2], [1, 1, 4]]
Expected: ['NO']
Received: ['YES']
```
여기서 Input이라고 돼있는 게 반례가 되는 거죠!

# 결론
이렇게 간단히 스트레스 테스팅과 그 예제를 알아봤어요.  
즐거운 디버깅 되시길 :)
