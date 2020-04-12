# library
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

with open("debug.txt") as input:
  ys = []
  xs = []
  rho = int(input.readline().strip())
  #next(input)
  for line in input:
    args = line.strip().split(",")
    xs.append(int(args[0]))
    ys.append(int(args[1]))
  xs = [math.log2(x) if x else 0 for x in xs]

  print("rho=" + str(rho) + " size=" + str(len(ys)))

  # stem function: first way
  #plt.stem(xs, ys)
  #plt.show()
  #plt.scatter(xs, ys, marker='o')
  #plt.show()


  fig = plt.figure()
  fig.subplots_adjust(bottom=0.2)
  ax = fig.add_subplot(111)

  plt.scatter(xs, ys, marker='o', s=0.2, color='violet', label = "AdaBoost")
  plt.show()

  #plt.scatter(xs,ys,alpha=0.1,cmap=cm.Paired)
  #plt.show()
  # stem function: If no X provided, a sequence of numbers is created by python:
  #plt.stem(values)
  #plt.show()

  # stem function: second way
  #(markerline, stemlines, baseline) = #plt.stem(xs, ys)
  #plt.setp(baseline, visible=False)
  #plt.show()

