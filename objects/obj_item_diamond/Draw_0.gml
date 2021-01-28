/// @description Insert description here
var _matrix_translate =        matrix_build(x,y, z,            0,0,0,                1,1,1);
var _matrix_rotate =            matrix_build(0,0,0,                                90,0,rotation,      1,1,1);
var _matrix_scale =                matrix_build(0,0,0,                                0,0,0,                16,16,16);
var _matrix_sr =                    matrix_multiply(_matrix_scale,_matrix_rotate);
var _matrix_final =                matrix_multiply(_matrix_sr,_matrix_translate);

matrix_set(matrix_world,_matrix_final);
vertex_submit(vb_diamond,pr_trianglelist,sprite_get_texture(spr_itemdiamond,0));
matrix_set(matrix_world,matrix_build_identity());