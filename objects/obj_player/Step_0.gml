/// @description movement and collision as well as sending packets

#region controls
//control initialization
var key_l = keyboard_check(vk_left);
var key_r = keyboard_check(vk_right);
var key_acc = keyboard_check(ord("Z"));
var key_sec = keyboard_check(ord("X"));
var key_drft = keyboard_check(ord("C"));
var key_back = keyboard_check(vk_down);

if keyboard_check(ord("R")) game_restart();
//changes turn speed based on your action
if action = 0 trnacc = 1.5;
if (action = 1 || action = 2) trnacc = 1.5;

//movement and attack/drifting
if key_drft && key_l && spd != 0 action = 1 else if key_drft && key_r && spd != 0 action = 2 else action = 0;
if key_l && action != 2 && spd != 0 angle += trnacc;
if key_r && action != 1 && spd != 0 angle -= trnacc;
if key_acc && !key_back spd += acc;
if key_acc && key_back spd -= acc;
#endregion

if spd > 0 spd -= dcc;
if spd < 0 spd += dcc;
if spd2 > 0 spd2 -= dcc;

if action = 0 && spd > mxspd spd = mxspd;
if action = 0 && spd < -mxspd spd = -mxspd;

if (action = 1 || action = 2) && spd > mxspddrft spd = mxspddrft;
if (action = 1 || action = 2) && spd < -mxspddrft spd = -mxspddrft;

#region player render angle
//changes player angle when you drift
if action = 0 && !key_l && !key_r r_angle = angle;
if action = 1 && spd != 0 r_angle = angle + 10;
if action = 2 && spd != 0 r_angle = angle - 10;

//when you turn changes angle
if action = 0 && spd != 0 && key_l && !key_r r_angle = angle + 5;
if action = 0 && spd != 0 && !key_l && key_r r_angle = angle - 5;
#endregion 

if (action = 1 || action = 2) chargeup = 1 else chargeup = 0;
if chargeup = 0 alarm[0] = 100;

if boost = 1 && action = 0 {
	spd2 = 10;
	boost_stop = 1;
}

if boost_stop = 0 alarm[1] = 70;

//idk it just like 
hsp = lengthdir_x(spd + spd2,angle);
vsp = lengthdir_y(spd + spd2,angle);
//zsp = lengthdir_x(,);

#region collision code
if place_meeting(x, y+vsp, par_solid)
{
while(!place_meeting(x, y+sign(vsp), par_solid)) //whilst the next pixel isn't a wall
{
y += sign(vsp);
}
vsp = 0;
}

if place_meeting(x+hsp, y, par_solid)
{
while(!place_meeting(x+sign(hsp), y, par_solid)) //whilst the next pixel isn't a wall
{
x += sign(hsp);
}
hsp = 0;
}

var slopeAngle = 60; //The threshhold where if the slope is steeper, the player will start to slide
var fast = false; //You’ll usually want to turn off the “fast” collision checking for important objects like the player
col = levelColmesh.displaceCapsule(x, y, z, 0, 0, 1, 48, 0, slopeAngle, fast);
if (col[6]) //If we’re touching ground
{
    x = col[0];
    y = col[1];
    z = col[2];
	ground = 1;
} else ground = 0;


  #endregion
var grv = -0.15;
if zsp < -5 zsp = -5;
if ground = 0 zsp -= grv;
if ground = 1 zsp = 0;
//applies movement
x += hsp;
y += vsp;
z += zsp;