l=5.;
Point(1) = {91.30000,00.00000,00.00000,l};
Point(2) = {121.30000,00.00000,00.00000,l};
Point(3) = {61.30000,00.00000,00.00000,l};
Point(4) = {91.30000,30.00000,00.00000,l};
Point(5) = {456.50000,00.00000,00.00000,l};
Point(6) = {486.50000,00.00000,00.00000,l};
Point(7) = {426.50000,00.00000,00.00000,l};
Point(8) = {456.50000,30.00000,00.00000,l};
Point(9) = {456.50000,316.27248,00.00000,l};
Point(10) = {486.50000,316.27248,00.00000,l};
Point(11) = {426.50000,316.27248,00.00000,l};
Point(12) = {456.50000,286.27248,00.00000,l};
Point(13) = {91.30000,316.27248,00.00000,l};
Point(14) = {121.30000,316.27248,00.00000,l};
Point(15) = {61.30000,316.27248,00.00000,l};
Point(16) = {91.30000,286.27248,00.00000,l};
Point(17) = {182.60000,158.13624,00.00000,l};
Point(18) = {212.60000,158.13624,00.00000,l};
Point(19) = {152.60000,158.13624,00.00000,l};
Point(20) = {182.60000,188.13624,00.00000,l};
Point(21) = {182.60000,128.13624,00.00000,l};
Point(22) = {365.20000,158.13624,00.00000,l};
Point(23) = {395.20000,158.13624,00.00000,l};
Point(24) = {335.20000,158.13624,00.00000,l};
Point(25) = {365.20000,188.13624,00.00000,l};
Point(26) = {365.20000,128.13624,00.00000,l};
Point(27) = {273.90000,00.00000,00.00000,l};
Point(28) = {328.90000,00.00000,00.00000,l};
Point(29) = {218.90000,00.00000,00.00000,l};
Point(30) = {273.90000,55.00000,00.00000,l};
Point(31) = {273.90000,316.27248,00.00000,l};
Point(32) = {328.90000,316.27248,00.00000,l};
Point(33) = {218.90000,316.27248,00.00000,l};
Point(34) = {273.90000,261.27248,00.00000,l};
Point(35) = {547.80000,158.13624,00.00000,l};
Point(36) = {492.80000,158.13624,00.00000,l};
Point(37) = {547.80000,213.13624,00.00000,l};
Point(38) = {547.80000,103.13624,00.00000,l};
Point(39) = {00.00000,158.13624,00.00000,l};
Point(40) = {55.00000,158.13624,00.00000,l};
Point(41) = {00.00000,213.13624,00.00000,l};
Point(42) = {00.00000,103.13624,00.00000,l};
Point(43) = {00.00000,00.00000,00.00000,l};
Point(44) = {547.80000,316.27248,00.00000,l};
Point(45) = {00.00000,316.27248,00.00000,l};
Point(46) = {547.80000,00.00000,00.00000,l};

Circle(1) = {15, 13, 16};
Circle(2) = {16, 13, 14};
Circle(3) = {33, 31, 34};
Circle(4) = {34, 31, 32};
Circle(5) = {11, 9, 12};
Circle(6) = {10, 9, 12};
Circle(7) = {37, 35, 36};
Circle(8) = {36, 35, 38};
Circle(9) = {6, 5, 8};
Circle(10) = {7, 5, 8};
Circle(11) = {29, 27, 30};
Circle(12) = {28, 27, 30};
Circle(13) = {3, 1, 4};
Circle(14) = {2, 1, 4};
Circle(15) = {42, 39, 40};
Circle(16) = {40, 39, 41};
Circle(17) = {19, 17, 20};
Circle(18) = {18, 17, 20};
Circle(19) = {19, 17, 21};
Circle(20) = {21, 17, 18};
Circle(21) = {24, 22, 25};
Circle(22) = {25, 22, 23};
Circle(23) = {23, 22, 26};
Circle(24) = {26, 22, 24};
Line(25) = {41, 45};
Line(26) = {45, 15};
Line(27) = {14, 33};
Line(28) = {32, 11};
Line(29) = {10, 44};
Line(30) = {44, 37};
Line(31) = {38, 46};
Line(32) = {46, 6};
Line(33) = {7, 28};
Line(34) = {29, 2};
Line(35) = {3, 43};
Line(36) = {43, 42};
Line Loop(37) = {25, 26, 1, 2, 27, 3, 4, 28, 5, -6, 29, 30, 7, 8, 31, 32, 9, -10, 33, 12, -11, 34, 14, -13, 35,
36, 15, 16};
Line Loop(38) = {17, -18, -20, -19};
Line Loop(39) = {21, 22, 23, 24};
Plane Surface(40) = {37, 38, 39};

