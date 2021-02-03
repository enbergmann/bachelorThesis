function l1NormOfJump4s = computeL1NormOfJump4s(currData, output)

  % extract necessary information from currData
  n4e = currData.n4e;
  n4s = currData.n4s;
  nrSides = currData.nrSides;
  s4e = currData.s4e;
  length4s = currData.length4s;
  e4s = currData.e4s;

  % extract necessary information from output
  u = output.u;

  % compute l1NormOfJump4s
  nodeValues4e = computeNodeValuesCR4e(s4e, u); 
  absNodeJumps4s = computeAbsNodeJumps4s(n4e, n4s, e4s, nrSides, nodeValues4e);

  l1NormOfJump4s = 1/4*length4s.*sum(absNodeJumps4s, 2);
    % for any edge F := conv({a, b}) = Tp \cap Tm, m := (a+b)/2,
    % u := u_{CR}, up := u|_Tp, and um := u|_Tm, and d := |up - um| it holds
    %   ||[u]_F||_{L^1(F)} = \int_F |up - um|(x) dx = \int_a^b d(x) dx 
    %                      = \int_a^m d(x) dx + \int_m^b d(x) dx
    % by continuity of the CR function u in m and affine linearity of up and
    % um the difference d = |up - um| satisfies d(m) = 0 and is affine linear
    % on both [a, m] and [m, b] with
    %   d(x) = (d(m) - d(a))/(m - a) (x - a) + d(a) 
    %        = -d(a)/(m - a) (x - a) + d(a)             
    % on [a, m] and
    %   d(x) = (d(b) - d(m))/(b - m) (x - m) + d(m)
    %        = d(b)/(b - m) (x - m)
    % on [m, b]
    % the midpoints rule, which is exact for affine linear functions, implies
    %   ||[u]_F||_{L^1(F)} = |m - a|d((m + a)/2) + |b - m|d((b + m)/2)
    % by definition of m it holds m - a = (b - a)/2 = b - m, hence
    %   ||[u]_F||_{L^1(F)} = |b - a|/2(d((m + a)/2) + d((b + m)/2))
    %                      = |b - a|/2(-d(a)/(m - a) ((m + a)/2 - a) + d(a) 
    %                                  + d(b)/(b - m) ((b + m)/2 - m))
    %                      = |b - a|/2(-d(a)/(m - a) (m - a)/2 + d(a) 
    %                                  + d(b)/(b - m) (b - m)/2)
    %                      = |b - a|/2(-d(a)/2 + d(a) + d(b)/2)
    %                      = |b - a|/2(d(a)/2 + d(b)/2)
    %                      = |b - a|/4(d(a) + d(b))
    % altogether
    %   ||[u]_F||_{L^1(F)} = |b - a|/4(|up(b) - um(b)| + |up(a) - um(a)|)
    %
    % and for any boundary edge F, by the midpoint rule and by u\in CR^1_0,
    % i.e. u(m) = 0
    %   ||[u]_F||_{L^1(F)} = \int_F |u|(x) dx
    %                      = \int_a^m |u|(x) dx + \int_m^b |u|(x) dx
    %                      = |m - a||u|((m + a)/2) + |b - m||u|((b + m)/2)
    %                      = |m - a|/2|u|(a) + |b - m|/2|u|(b)
    %                      = |b - a|/2(|u|(a) + |u|(b))
    %                      = |b - a||u|(a)                     (|u|(a)=|u|(b))
    %   TODO ask if this is correct and if boundary sides work this way
    %        i.e. are not just 0
    %        ALSO check if this is implemented (probably not the case)
    %        ALSO is this really correct for boundary edges in estimate?
end
