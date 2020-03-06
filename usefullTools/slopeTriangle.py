from math import log10

print("\nLet m = k/l = (lg(y2)-lg(y1))/(lg(x2)-lg(x1)).\n")
k = float(input("k = "))
l = float(input("l = "))
m = k/l
coord =\
  input("\nPlease enter the coordinate you want to compute (y2, y1, x2, x1): ")

coordFlag = True

if coord == "y2":
  y1 = float(input("\ny1 = "))
  x2 = float(input("x2 = "))
  x1 = float(input("x1 = "))
  resLg = m*log10(x2/x1) + log10(y1)
elif coord == "y1":
  y2 = float(input("\ny2 = "))
  x2 = float(input("x2 = "))
  x1 = float(input("x1 = "))
  resLg = -m*log10(x2/x1) + log10(y2)
elif coord == "x2":
  y2 = input("\ny2 = ")
  y1 = input("y1 = ")
  x1 = input("x1 = ")
  resLg = log10(y2/y1)/m + log10(x1)
elif coord == "x1":
  y2 = input("\ny2 = ")
  y1 = input("y1 = ")
  x2 = input("x2 = ")
  resLg = -log10(y2/y1)/m + log10(x2)
else:
  coordFlag = False  
  print("\nFlawed input: ", coord, "\n")
  
if coordFlag:
  print("\n", coord, " = ", 10**resLg, "\n")
