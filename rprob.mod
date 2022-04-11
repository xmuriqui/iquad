#AMPL model file generated by minlp problem library, developed by Wendel Melo (Ph. D. student, Federal University of Rio de Janeiro, Brazil)


var x0;
var x1;
var x2;
var x3;
var x4 >= 0 <= 1.4142135623731;
var x5 >= 0 <= 1.4142135623731;
var x6 >= 0 <= 1.4142135623731;
var x7 >= 0 <= 1.4142135623731;
var x8 >= 0 <= 1.4142135623731;
var x9 >= 0 <= 1.4142135623731;
var x10 >= 0 <= 1.4142135623731;
var x11 >= 0 <= 1.4142135623731;
var x12 integer >= 0 <= 1;
var x13 integer >= 0 <= 1;
var x14 integer >= 0 <= 1;
var x15 integer >= 0 <= 1;
var x16 integer >= 0 <= 1;
var x17 integer >= 0 <= 1;
var x18 integer >= 0 <= 1;
var x19 integer >= 0 <= 1;
var x20 >= 0 <= 1.4142135623731;
var x21;
var x22;
var x23;
var x24;
var x25;
var x26;
var x27;
var x28;

minimize obj: +1*x20 +0.5*x21 +0.5*x22 +0.5*x23 +0.5*x24 +0.5*x25 +0.5*x26 +0.5*x27 +0.5*x28 +0.2500005*x4*x4 +0.2500005*x5*x5 +0.2500005*x6*x6 +0.2500005*x7*x7 +0.2500005*x8*x8 +0.2500005*x9*x9 +0.2500005*x10*x10 +0.2500005*x11*x11 +0.5*x12*x4 +0.2500005*x12*x12 +0.5*x13*x5 +0.2500005*x13*x13 +0.5*x14*x6 +0.2500005*x14*x14 +0.5*x15*x7 +0.2500005*x15*x15 +0.5*x16*x8 +0.2500005*x16*x16 +0.5*x17*x9 +0.2500005*x17*x17
 +0.5*x18*x10 +0.2500005*x18*x18 +0.5*x19*x11 +0.2500005*x19*x19;

#constraints 0 1 2 3 4 5 6 7 8 have nonlinear components, but, unfortunatelly, those nonlinear components cannot be represented in this file.
subject to cons0:  0 <= -0;
subject to cons1:  0 <= -0;
subject to cons2:  0 <= -0;
subject to cons3:  0 <= -0;
subject to cons4:  0 <= -0;
subject to cons5:  0 <= -0;
subject to cons6:  0 <= -0;
subject to cons7:  0 <= -0;
subject to cons8:  0 <= -0;
subject to cons9: 1 <=  +1*x12 +1*x13 <= 1;
subject to cons10: 1 <=  +1*x14 +1*x15 <= 1;
subject to cons11: 1 <=  +1*x16 +1*x17 <= 1;
subject to cons12: 1 <=  +1*x18 +1*x19 <= 1;
subject to cons13: 2 <=  +1*x12 +1*x14 +1*x16 +1*x18 <= 2;
subject to cons14: 2 <=  +1*x13 +1*x15 +1*x17 +1*x19 <= 2;
subject to cons15: -0.999999999999072 <=  -0.707106781186547*x8 +0.707106781186547*x16 <= 0.70003568289672;
subject to cons16: -0.999999999932297 <=  -0.707106781186547*x10 +0.707106781186547*x18 <= 0.700035669197126;
subject to cons17: -0.700035669179761 <=  +0.707106781186547*x4 -0.707106781186547*x12 <= 0.999999999983433;
subject to cons18: -0.70003568300068 <=  +0.707106781186547*x6 -0.707106781186547*x14 <= 0.999999999999071;
subject to cons19: -0.700035683004388 <=  +0.707106781186547*x9 -0.707106781186547*x17 <= 0.999999999999071;
subject to cons20: -0.700035669197172 <=  +0.707106781186547*x11 -0.707106781186547*x19 <= 0.999999999932732;
subject to cons21: -0.999999999983433 <=  -0.707106781186547*x5 +0.707106781186547*x13 <= 0.70003566917976;
subject to cons22: -0.999999999999062 <=  -0.707106781186547*x7 +0.707106781186547*x15 <= 0.70003568301719;