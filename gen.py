#!/usr/bin/python3

import sys
import random

LB = -100000000
UB =  100000000

try:
  n = int(sys.argv[1])
except:
  print("Usage: ./gen.py <side length>")
  quit()

with open("input/{}x{}.txt".format(str(n), str(n)), "w") as f:
  f.write(str(n) + '\n')
  
  nums = random.choices(range(LB, UB), k=n*n)

  for row in range(n):
    for col in range(n):
      f.write(str(nums.pop()) + ' ')
    f.write('\n')
