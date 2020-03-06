from math import log10

print("\nLet m = k/l = (lg(y)-lg(yn))/(lg(x)-lg(xn)).\n")
k = float(input("k = "))
l = float(input("l = "))
m = k/l
coord =\
  input("\nPlease enter the coordinate you want to compute (yn, xn): ")

coordFlag = True

if coord == "yn":
  x2 = float(input("\nx2 = "))
  y2 = float(input("y2 = "))
  xn = float(input("xn = "))
  resLg = -m*log10(x2/xn) + log10(y2)
elif coord == "xn":
  x2 = float(input("\nx2 = "))
  y2 = float(input("y2 = "))
  yn = float(input("yn = "))
  resLg = -log10(y2/yn)/m + log10(x2)
else:
  coordFlag = False  
  print("\nFlawed input: ", coord, "\n")
  
if coordFlag:
  print("\n", coord, " = ", 10**resLg, "\n")
