l=1.0;
Point(1) = {0,0,0,l};
Point(2) = {1,0,0,l};
Point(3) = {1,1,0,l};
Point(4) = {0,1,0,l};
Point(5) = {0.25,0.25,0,l};
Point(6) = {0.75,0.25,0,l};
Point(7) = {0.75,0.75,0,l};
Point(8) = {0.25,0.75,0,l};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line(5) = {8, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(10) = {9};