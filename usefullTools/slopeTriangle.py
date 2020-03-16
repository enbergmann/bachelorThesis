from math import log10

print("\nLet m = k/l = (lg(y)-lg(yn))/(lg(x)-lg(xn)).\n")
k = float(input("k = "))
l = float(input("l = "))
m = k/l
coord =\
  input("\nPlease enter the coordinate you want to compute (yn, xn): ")

coordFlag = True

if coord == "yn":
  x = float(input("\nx = "))
  y = float(input("y = "))
  xn = float(input("xn = "))
  resLg = -m*log10(x/xn) + log10(y)
elif coord == "xn":
  x = float(input("\nx = "))
  y = float(input("y = "))
  yn = float(input("yn = "))
  resLg = -log10(y/yn)/m + log10(x)
else:
  coordFlag = False  
  print("\nFlawed input: ", coord, "\n")
  
if coordFlag:
  print("\n", coord, " = ", 10**resLg, "\n")
