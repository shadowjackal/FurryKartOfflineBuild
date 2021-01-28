/// @description 
var camera = camera_get_active();
var camera_distance = 192;


#region most variables
var xto = obj_player.x;
var yto = obj_player.y;
var zto = obj_player.z;
var xfrom = xto - camera_distance * dcos(obj_player.angle);
var yfrom = yto + camera_distance * dsin(obj_player.angle);
var zfrom = zto + camera_distance * dsin(-20);
#endregion

#region camera direction and such
camera_set_view_mat(camera, matrix_build_lookat(xfrom, yfrom, zfrom, xto, yto, zto, 0, 0, 1));
camera_set_proj_mat(camera, matrix_build_projection_perspective_fov(60, window_get_width() / window_get_height(), 1, 32000));
camera_apply(camera);
#endregion 

#region World Matrix Build wall
//Order matters because of the way translations are handled (centered on 0,0)
var _matrix_translate =        matrix_build(obj_player.x,obj_player.y, obj_player.z,            0,0,0,                1,1,1);
var _matrix_rotate =            matrix_build(0,0,0,                                -obj_player.z_angle + 90,0,obj_player.r_angle-90,      1,1,1);
var _matrix_scale =                matrix_build(0,0,0,                                0,0,0,                16,16,16);
var _matrix_sr =                    matrix_multiply(_matrix_scale,_matrix_rotate);
var _matrix_final =                matrix_multiply(_matrix_sr,_matrix_translate);

matrix_set(matrix_world,_matrix_final);
vertex_submit(vb_player,pr_trianglelist,sprite_get_texture(spr_player_texture,0));
matrix_set(matrix_world,matrix_build_identity());
#endregion

#region World Matrix Build wall
//Order matters because of the way translations are handled (centered on 0,0)
var _matrix_translate =        matrix_build(0,0, 0,            0,0,0,                1,1,1);
var _matrix_rotate =            matrix_build(0,0,0,                                 0,0,90,      1,1,1);
var _matrix_scale =                matrix_build(0,0,0,                                0,0,0,                1,1,1);
var _matrix_sr =                    matrix_multiply(_matrix_scale,_matrix_rotate);
var _matrix_final =                matrix_multiply(_matrix_sr,_matrix_translate);

matrix_set(matrix_world,_matrix_final);
vertex_submit(vb_level,pr_trianglelist,-1);
matrix_set(matrix_world,matrix_build_identity());
#endregion

