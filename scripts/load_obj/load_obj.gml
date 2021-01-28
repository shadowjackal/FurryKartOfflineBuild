///@func load_obj_file(obj_file_name,[mtl_file_name])
///@desc Loads an object file
///@arg obj_file_name
///@arg [mtl_file_name]
function load_obj() {
	var _obj_file_name = argument[0];
	if(argument_count > 1) { var _mtl_file_name = argument[1]; }

#region Setup
	var _timer = get_timer();
	var _obj_file = file_text_open_read(_obj_file_name);
	if(argument_count > 1) { var _mtl_file = file_text_open_read(_mtl_file_name); }

	vertex_format_begin();
	vertex_format_add_position_3d();
	vertex_format_add_normal();
	vertex_format_add_texcoord();
	vertex_format_add_color();
	var _vformat = vertex_format_end();
	var _model = vertex_create_buffer();
	vertex_begin(_model,_vformat);

	var _vertex_x = ds_list_create();
	var _vertex_y = ds_list_create();
	var _vertex_z = ds_list_create();
	var _vertex_nx = ds_list_create();
	var _vertex_ny = ds_list_create();
	var _vertex_nz = ds_list_create();
	var _vertex_u = ds_list_create();
	var _vertex_v = ds_list_create();

	var _mtl_alpha = ds_map_create();
	var _mtl_color = ds_map_create();
	ds_map_set(_mtl_color,"none",c_white);
	ds_map_set(_mtl_alpha,"none",1);
	var _active_mtl = "none";
#endregion

#region Mtl File Parsing
	if(argument_count > 1) {
		while(!file_text_eof(_mtl_file)) {
			var _line = file_text_readln(_mtl_file);
			var _index = 0;
			var _terms = [""];
			_terms[string_count(_line," ")] = "";
			for(var i=1; i<=string_length(_line); i++) {
				if(string_char_at(_line,i) == " ") {
					_index++;
					_terms[_index] = "";
				} else {
					if(!any(string_char_at(_line,i),"\r","\n")) { _terms[_index] += string_char_at(_line,i); }
				}
			}
	
			var _color = c_white;
			var _alpha = 1;
			switch(_terms[0]) {
				case "newmtl":
					_name = _terms[1];
				break;
				case "Kd":
					var _red = real(_terms[1]) * 255;
					var _green = real(_terms[2]) * 255;
					var _blue = real(_terms[3]) * 255;
					_color = make_color_rgb(_red,_green,_blue);
					ds_map_set(_mtl_color,_name,_color);
				break;
				case "d":
					_alpha = real(_terms[1]);
					ds_map_set(_mtl_alpha,_name,_alpha);
				break;
				default:
					//show_debug_message("WARNING: Attempting to access an unsupported mtl type in mtl file " + _mtl_file_name);
				break;
			}
		}
	}
#endregion

#region Obj File Parsing
	while(!file_text_eof(_obj_file)) {
		var _line = file_text_read_string(_obj_file);
		file_text_readln(_obj_file);
		var _index = 0;
		var _terms = [""];
		_terms[string_count(_line," ")] = ""; //Presize our array, a very slight optimization
		for(var i=1; i<=string_length(_line); i++) {
			if(string_char_at(_line,i) == " ") {
				_index++;
				_terms[_index] = "";
			} else {
				if(!any(string_char_at(_line,i),"\r","\n")) { _terms[_index] += string_char_at(_line,i);	}
			}
		}
	
		switch(_terms[0]) {
			case "v":
				ds_list_add(_vertex_x,real(_terms[1]));
				ds_list_add(_vertex_y,real(_terms[2]));
				ds_list_add(_vertex_z,real(_terms[3]));
			break;
			case "vt":
				ds_list_add(_vertex_u,real(_terms[1]));
				ds_list_add(_vertex_v,real(_terms[2]));
			break;
			case "vn":
				ds_list_add(_vertex_nx,real(_terms[1]));
				ds_list_add(_vertex_ny,real(_terms[2]));
				ds_list_add(_vertex_nz,real(_terms[3]));
			break;
			case "f":
				for(var k=1; k<=3; k++) {
					var _index = 0;
					var _data = [""];
					_data[string_count(_terms[k],"/")] = "";
					for(var a=1; a<=string_length(_terms[k]); a++) {
						if(string_char_at(_terms[k],a) == "/") {
							_index++;
							_data[_index] = "";
						} else {
							_data[_index] += string_char_at(_terms[k],a);	
						}
					}
				
					var xx = ds_list_find_value(_vertex_x,real(_data[0]) - 1);
					var yy = ds_list_find_value(_vertex_y,real(_data[0]) - 1);
					var zz = ds_list_find_value(_vertex_z,real(_data[0]) - 1);
					var u  = ds_list_find_value(_vertex_u,real(_data[1]) - 1);
					var v = 1 - ds_list_find_value(_vertex_v,real(_data[1]) - 1);
					var nx = ds_list_find_value(_vertex_nx,real(_data[2]) - 1);
					var ny = ds_list_find_value(_vertex_ny,real(_data[2]) - 1);
					var nz = ds_list_find_value(_vertex_nz,real(_data[2]) - 1);
				
					var _color = c_white;
					var _alpha = 1;
					if(ds_map_exists(_mtl_color,_active_mtl)) { _color = _mtl_color[? _active_mtl]; }
					if(ds_map_exists(_mtl_alpha,_active_mtl)) { _alpha = _mtl_alpha[? _active_mtl]; }
				
					vertex_position_3d(_model,xx,yy,zz);
					vertex_normal(_model,nx,ny,nz);
					vertex_texcoord(_model,u,v);
					vertex_color(_model,_color,_alpha);
				}
			break;
			case "usemtl":
				_active_mtl = _terms[1];
			break;
			default:
				//show_debug_message("WARNING: Attempting to access an unsupported obj type in obj file " + _obj_file_name);
			break;
		}
	}
	vertex_end(_model);
#endregion

#region Cleanup
	ds_list_destroy(_vertex_x);
	ds_list_destroy(_vertex_y);
	ds_list_destroy(_vertex_z);
	ds_list_destroy(_vertex_nx);
	ds_list_destroy(_vertex_ny);
	ds_list_destroy(_vertex_nz);
	ds_list_destroy(_vertex_u);
	ds_list_destroy(_vertex_v);

	ds_map_destroy(_mtl_alpha);
	ds_map_destroy(_mtl_color);

	file_text_close(_obj_file);
	if(argument_count > 1) { file_text_close(_mtl_file); }

	var _end_timer = (get_timer() - _timer) / 1000;
	show_debug_message("Load time for " + _obj_file_name + " is: " + string(_end_timer) + " ms");
#endregion


	return _model;


}
