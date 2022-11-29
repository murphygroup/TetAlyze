function output_area = calculate_triangle_area(x1,y1,z1,x2,y2,z2,x3,y3,z3)
%CALCULATE_TRIANGLE_AREA Summary of this function goes here
%   Detailed explanation goes here
    a=((x1-x2)^2.0 + (y1-y2)^2.0 + (z1-z2)^2.0)^0.5;
    b=((x3-x2)^2.0 + (y3-y2)^2.0 + (z3-z2)^2.0)^0.5;
    c=((x3-x1)^2.0 + (y3-y1)^2.0 + (z3-z1)^2.0)^0.5;
    s = (a + b + c) / 2   
    output_area = (s*(s-a) * (s-b)*(s-c))^0.5        

    

end

