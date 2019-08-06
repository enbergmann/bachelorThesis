function [ Phi, GradPhi ,nCP, nIP , c4Base ] = getP3ShapeFunctions( ) 
%%GETP3SHAPEFUNCTION Create P3 basis functions and their gradients.
%   GETP3SHAPEFUNCTION() returns a structure array containing the
%   function handle of the P3 basis functions on the reference triangle.
%
%   See also TRANSFORMAFFINE.
%
% output    Phi         - basis on the reference triangle for total degree 
%                         3; order: base functions in vertices,on first side,
%                         second side, third side, of inner points
%           GradPhi     - gradients of Phi (in the same order)
%           nCP         - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIP         - number of base functions in inner points

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a


phi3 = @(x) ( 1/2 * (1 - x(:,1) - x(:,2)) .* (3 * (1 - x(:,1) - x(:,2)) -1 ) .* ( 3 * (1 - x(:,1) - x(:,2))-2 ) ) ;%B2
phi1 = @(x) ( 1/2 * x(:,1) .* (3 * x(:,1) -1 )                .* ( 3 * x(:,1) -2) );%B3
phi2 = @(x) ( 1/2 * x(:,2) .* (3 * x(:,2) -1 )                .* ( 3 * x(:,2) -2) );%B4

phi8 = @(x) (9/2 *  x(:,1)               .* (3 *  x(:,1) - 1)               .*  x(:,2) );%B8
phi9 = @(x) (9/2 *  x(:,2)               .* (3 *  x(:,2) - 1)               .*  x(:,1) );%B10

phi4 = @(x) ( 9/2 *  x(:,2)               .* (3 *  x(:,2) - 1)               .* (1 - x(:,1) - x(:,2)) );%B9
phi5 = @(x) ( 9/2 * (1 - x(:,1) - x(:,2)) .* (3 * (1 - x(:,1) - x(:,2)) - 1) .*  x(:,2) );%B6

phi6 = @(x) ( 9/2 * (1 - x(:,1) - x(:,2)) .* (3 * (1 - x(:,1) - x(:,2)) - 1) .*  x(:,1) );%B5
phi7 = @(x) ( 9/2 *  x(:,1)               .* (3 *  x(:,1) - 1)               .* (1 - x(:,1) - x(:,2)) );%B7

phi10 = @(x) ( 27* (1 - x(:,1) - x(:,2)) .*  x(:,1) .*  x(:,2)); %B1





gradPhi3 = @(x) [- 1/2 * (-2 + 3 * (1 -x(:,1)-x(:,2))).* (-1 + 3 * (1 -x(:,1)-x(:,2))) - 3/2 * (-2 + 3 * (1 -x(:,1)-x(:,2))) .* (1 -x(:,1)-x(:,2)) - 3/2 * (-1 + 3 * (1 -x(:,1)-x(:,2))) .* (1 -x(:,1)-x(:,2)) ,...
                 - 1/2 * (-2 + 3 * (1 -x(:,1)-x(:,2))).* (-1 + 3 * (1 -x(:,1)-x(:,2))) - 3/2 * (-2 + 3 * (1 -x(:,1)-x(:,2))) .* (1 -x(:,1)-x(:,2)) - 3/2 * (-1 + 3 * (1 -x(:,1)-x(:,2))) .* (1 -x(:,1)-x(:,2)) ];
            %checked
gradPhi1 = @(x) [ 3/2 * x(:,1) .* (-2 + 3 * x(:,1))...
                + 3/2 * x(:,1) .* (-1 + 3 * x(:,1)) + 1/2 * (-2 + 3 * x(:,1)) .* (-1 + 3 * x(:,1)), zeros(size(x,1),1)]; 
            %checked
gradPhi2 = @(x) [ zeros(size(x,1),1), 3/2 * x(:,2) .* (-2 + 3 * x(:,2)) + 3/2 * x(:,2) .* (-1 + 3 * x(:,2)) + 1/2 * (-2 + 3 * x(:,2)).* (-1 + 3 * x(:,2)) ]; 
            %checked
gradPhi8 = @(x) [ 27/2 * x(:,1) .* x(:,2) + 9/2 * (-1 + 3* x(:,1)) .* x(:,2) ,  9/2 * x(:,1) .* (-1 + 3 * x(:,1))  ]; 
            %checked

gradPhi9 = @(x) [ 9/2 * x(:,2) .* (-1 + 3 * x(:,2)), 27/2 * x(:,1) .* x(:,2) + 9/2 * x(:,1) .* (-1 + 3 * x(:,2))]; 
            %checked

gradPhi4 = @(x) [-9/2 * x(:,2) .* (-1 + 3 * x(:,2)) , -9/2 * x(:,2) .* (-1 + 3 * x(:,2)) + 27/2 * x(:,2) .* (1 - x(:,1) - x(:,2)) + 9/2 * (-1 + 3 * x(:,2)) .* (1- x(:,1) - x(:,2))  ]; 
            %checked

gradPhi5 = @(x) [-9/2 * x(:,2) .* (-1 + 3 * (1 - x(:,1) - x(:,2))) - 27/2 * x(:,2) .* (1 - x(:,1) - x(:,2)) ,...
                 -9/2 * x(:,2) .* (-1 + 3 * (1 - x(:,1) - x(:,2))) - 27/2 * x(:,2) .* (1 - x(:,1) - x(:,2)) + 9/2 * (-1 + 3 * (1- x(:,1) - x(:,2))) .* (1- x(:,1) - x(:,2))  ]; 
             %checked            
gradPhi6 = @(x) [-9/2 * x(:,1) .* (-1 + 3 * (1 - x(:,1) - x(:,2))) ...
                 - 27/2 *x(:,1) .* (1 - x(:,1) - x(:,2)) + 9/2 * (-1 + 3 * (1- x(:,1) - x(:,2))) .* (1- x(:,1) - x(:,2)), -9/2 * x(:,1) .* (-1 + 3 * (1 - x(:,1) - x(:,2))) - 27/2 * x(:,1) .* (1 - x(:,1) - x(:,2))] ;

gradPhi7 = @(x) [-9/2 * x(:,1) .* (-1 + 3 * x(:,1)) + 27/2 * x(:,1) .* (1 - x(:,1) - x(:,2)) + 9/2 * (-1 + 3 * x(:,1)) .* (1- x(:,1) - x(:,2)) , - 9/2 * x(:,1) .* (-1 + 3 * x(:,1))  ];


gradPhi10 = @(x)[ - 27 * x(:,1) .* x(:,2) + 27 * (1 - x(:,1) - x(:,2)) .* x(:,2), ...
                   27 * x(:,1) .* (1 - x(:,1) - x(:,2)) - 27 * x(:,1) .* x(:,2)  ];

Phi = {phi1, phi2, phi3, phi4, phi5, phi6,phi7,phi8,phi9,phi10};
GradPhi = {gradPhi1, gradPhi2, gradPhi3, gradPhi4, gradPhi5, gradPhi6, ...
    gradPhi7, gradPhi8, gradPhi9, gradPhi10};

nCP = 1; nIP = 1; c4Base = [1/3, 2/3]; 
end

