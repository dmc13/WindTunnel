basin_x = 2.2;
basin_y = 0.41;
site_x = 0.01;
site_y = 0.01;
rad = 0.05;
n = 50;
element_size = 0.005;
element_size_coarse = 0.05;
circle_x = 0.2;

Point(1) = {0, 0, 0, element_size};
Point(2) = {basin_x, 0, 0, element_size_coarse};
Point(3) = {0, basin_y, 0, element_size};
Point(4) = {basin_x, basin_y, 0, element_size_coarse};

//Point(5) = {circle_x, (basin_y/2-rad)/2-site_y/2, 0, element_size_coarse};
//Extrude{site_x, 0, 0} { Point{5}; Layers{site_x/element_size_coarse}; }
//Extrude{0, site_y, 0} { Line{1}; Layers{site_y/element_size_coarse}; }

Line(6) = {1, 2};
Line(7) = {2, 4};
Line(8) = {4, 3};
Line(9) = {3, 1};
//Line Loop(10) = {9, 6, 7, 8};
//Line Loop(11) = {3, 2, -4, -1};

Point(9) = {circle_x, basin_y/2, 0, element_size};
Point(10) = {circle_x, basin_y/2-rad, 0, element_size};
Point(11) = {circle_x, basin_y/2+rad, 0, element_size};
Circle(14) = {11, 9, 10};
Circle(15) = {10, 9, 11};
//Line Loop(16) = {14, 15};
//Plane Surface(17) = {10, 11, 16};

Physical Line(2) = {7};
Physical Line(1) = {9};
Physical Line(3) = {8, 6, 14, 15};
//Physical Surface(18) = {17};
Line Loop(16) = {8, 9, 6, 7};
Line Loop(17) = {15, 14};
Plane Surface(18) = {16, 17};
Physical Surface(18) = {18};
