#!/usr/bin/python3

import sys
import random

LB = -1000000000
UB = 10000000000

try:
  n = int(sys.argv[1])
except:
  print("Usage: ./gen.py <side length>")
  quit()

with open("input/{}x{}.txt".format(str(n), str(n)), "w") as f:
  f.write(str(n) + '\n')
  
  nums = random.sample(range(LB, UB), k=n*n)

  for row in range(n):
    for col in range(n):
      f.write(str(nums.pop()) + ' ')
    f.write('\n')
