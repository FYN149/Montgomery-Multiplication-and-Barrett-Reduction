import random
from Crypto.Util.number import getPrime, inverse

# 蒙哥马利约简
def montgomery_reduce(a, n, n_inv, r):
    # 计算中间变量q
    q = ((a & (r - 1)) * n_inv) & (r - 1)
    # 计算约简后的a
    a = (a + q * n) >> (r.bit_length() - 1)
    # 如果a大于等于模数n，则减去n
    if a >= n:
        a -= n
    return a

# 蒙哥马利乘法
def montgomery_multiply(a_bar, b_bar, n, n_inv, r):
    # 计算乘积
    t = a_bar * b_bar
    # 对乘积进行蒙哥马利约简
    return montgomery_reduce(t, n, n_inv, r)

# 巴雷特约简
def barrett_reduction(x, n):
    # 计算模数n的位数
    k = n.bit_length()
    # 计算2的2k次幂
    r = 1 << (k * 2)
    # 计算n的逆元（使用Crypto库的inverse函数）
    n_inv = inverse(n, r)
    # 计算中间变量q2
    q2 = ((x * n_inv) >> k) >> k
    # 计算约简后的结果
    r = x - q2 * n
    # 如果结果大于等于n，则减去n
    if r >= n:
        r -= n
    return r

# 生成两个随机的2048位数和一个2048位的素数模数
a = random.randint(2**2047, 2**2048 - 1)
b = random.randint(2**2047, 2**2048 - 1)
n = getPrime(2048)

# 初始化蒙哥马利乘法的参数
r = 1 << n.bit_length()  # r是大于n的最小2的幂
n_inv = inverse(n, r)  # n的模逆
# 将a和b转换为蒙哥马利形式
a_bar = (a * r) % n
b_bar = (b * r) % n

# 执行蒙哥马利乘法
montgomery_result = montgomery_multiply(a_bar, b_bar, n, n_inv, r)
# 将结果转换回普通形式
# 但为了与直接计算的结果对比，我们可以再次约简以确保一致性（实际上这一步是多余的）
montgomery_result_normal = montgomery_reduce(montgomery_result * n_inv, n, n_inv, r) % n

# 执行巴雷特约简
# 为了对比，我们直接对a*b的结果进行巴雷特约简
barrett_result = barrett_reduction(a * b, n)

# 输出结果
print("Answer")
print((a * b) % n)
print("Montgomery Multiplication (Normal Form)")
print(montgomery_result_normal)
print("Barrett Reduction ")
print(barrett_result)