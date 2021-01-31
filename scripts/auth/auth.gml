// Script assets have changed for v2.3.0 see
// https://help.yoyogames.com/hc/en-us/articles/360005277377 for more information
function auth(socket, buffer){
	buffer_seek(buffer, buffer_seek_start, 0);
	
	// ID 0x01
	buffer_write(buffer, buffer_u8, 0x01);
	
	// Username
	buffer_write(buffer, buffer_string, obj_player.username);
	
	network_send_packet(socket, buffer, buffer_tell(buffer));
}