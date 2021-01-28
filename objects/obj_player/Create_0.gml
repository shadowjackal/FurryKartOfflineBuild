/// @description Insert description here
// You can write your code in this editor
acc = 0.5;
dcc = 0.125;
trnacc = 1.5;
spd = 0;
spd2 = 0;
zsp = 0;
angle = 0;
r_angle = 0;
mxspd = 30;
mxspddrft = 15;
action = 0;
z = -10;
hsp = 0;
vsp = 0;
boost = 0;
boost_stop = 0;
item = 0;
z_angle = 0;
ground = 0;
rtate_dcc = 0.5;
rtatelimit = 12;
rtatelimitdrift = 20;
plusratate = 0;

var _matrix_translate =        matrix_build(0,0, 0,            0,0,0,                1,1,1);
var _matrix_rotate =            matrix_build(0,0,0,                                 0,0,90,      1,1,1);
var _matrix_scale =                matrix_build(0,0,0,                                0,0,0,                1,1,1);
var _matrix_sr =                    matrix_multiply(_matrix_scale,_matrix_rotate);
var _matrix_final = matrix_multiply(_matrix_sr,_matrix_translate);
var M = matrix_multiply(matrix_build(0, 0, 0, 90, 0, 0, 1, 1, -1), _matrix_final);

levelColmesh = new colmesh();
levelColmesh.addMesh("lebeltest.obj", M);
var regionSize = 1000; //120 is a magic number I chose that fit well for my player size and level complexity. It may have to be different for your game!
levelColmesh.subdivide(regionSize);