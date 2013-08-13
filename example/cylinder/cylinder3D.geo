
b=21.7000000;
c=12.0000000;
l=2.00000;
Point(1) = {-b, -b, 0, l};
Point(2) = {b, -b, 0, l};
Point(3) = {b, b, 0, l};
Point(4) = {-b, b, 0, l};
Point(5) = {0.0, 0.0, 0, l};
Point(6) = {-c, 0.0, 0, l};
Point(7) = {c, 0.0, 0, l};
Point(8) = {0.0, -c, 0, l};
Point(9) = {0.0, c, 0, l};

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
Extrude {0, 0, 2*b} {
  Surface{12};
}
