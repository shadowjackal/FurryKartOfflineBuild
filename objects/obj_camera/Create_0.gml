/// @description vertex buffers and other useful stuff
audio_stop_all();
//audio_play_sound(FurryKart_Racing_Ver_3,1,true);


//this area is for the set up for 3d camera and drawing objects
#region gpu_set_z
gpu_set_ztestenable(true);
gpu_set_zwriteenable(true);
#endregion


#region vertex format and set up
vertex_format_begin();
vertex_format_add_position_3d();
vertex_format_add_normal();
vertex_format_add_texcoord();
vertex_format_add_color();
vertex_format = vertex_format_end();

vbuffer = vertex_create_buffer();
vertex_begin(vbuffer, vertex_format); 
#endregion

vb_player = vertex_create_buffer();
vb_player = load_obj("furrykart.obj","furrykart.mtl");

vb_level = vertex_create_buffer();
vb_level = load_obj("lebeltest.obj","lebeltest.mtl");


