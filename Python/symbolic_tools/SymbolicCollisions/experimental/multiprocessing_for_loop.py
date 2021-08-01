import time
# from multiprocessing import Pool
from multiprocessing.pool import ThreadPool as Pool
import numpy as np

first_it = [1, 2, 3, 4]
second_it = [5, 6, 7, 8]

serial_result = list(map(pow, first_it, second_it))
print(serial_result)

p = Pool(2)
# parallel_result = p.map(np.sqrt, first_it)
# p.close()
# p.join()

# is there a variant of pool.map which support multiple arguments?
# Python 3.3 includes pool.starmap() method:
parallel_result = p.starmap(pow, zip(first_it, second_it))
p.close()
p.join()
print(parallel_result)


N = 2


# for i in range(N):
#     for j in range(N):
#         for k in range(N):
#             value = 100*i+10*j+k
#             T[i][j][k] = value
# x_array = np.arange(0, N, 1)
T = np.zeros((N, N, N))

def fun(i, j, k):
    value = 100*i+10*j+k
    return value
    # T[i][j][k] = value


p = Pool(2)

parallel_result = p.starmap(fun, zip(range(N), range(N), range(N)))
p.close()
p.join()
print(parallel_result)
