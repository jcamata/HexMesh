lc = 0.5;
Point(1) = {20.0 , 0.0, 2.0, lc};
Point(2) = {0.5  , 0.0, 2.0, lc};
Point(3) = {-0.5 , 0.0, -2.0,lc};
Point(4) = {-20.0, 0.0, -2.0, lc};


BSpline(1) = {1, 2, 3, 4};


Extrude {0, 20, 0} {
  Line{1};
}

Line Loop(6) = {2, -4, -1, 3};
Ruled Surface(7) = {6};
Coherence;
