l=0.05;
Point(1) = {-1., -1., 0, l};
Point(2) = {1., -1., 0, l};
Point(3) = {1., 1., 0, l};
Point(4) = {-1., 1., 0, l};
Point(5) = {0.0, 0.0, 0, l};
Point(6) = {-0.5, 0.0, 0, l};
Point(7) = {0.5, 0.0, 0, l};
Point(8) = {0.0, -0.5, 0, l};
Point(9) = {0.0, 0.5, 0, l};

Line(1) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line(5) = {3, 4};
Circle(6) = {6, 5, 8};
Circle(7) = {8, 5, 7};
Circle(8) = {7, 5, 9};
Circle(9) = {9, 5, 6};
Line Loop(10) = {5, 1, 3, 4};
Line Loop(11) = {8, 9, 6, 7};
Plane Surface(12) = {10, 11};
