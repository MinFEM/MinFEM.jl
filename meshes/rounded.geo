//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0.1};
//+
Physical Surface("domain", 10001) = {1};
//+
Physical Line("bottom", 1001) = {1};
//+
Physical Line("right", 1002) = {3};
//+
Physical Line("top", 1003) = {5};
//+
Physical Line("left", 1004) = {7};
//+
