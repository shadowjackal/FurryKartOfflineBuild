online_mode = 0;

address = "127.0.0.1";
port = 6621;

if (online_mode) {
	socket = network_create_socket(network_socket_tcp);
	network_connect(socket, address, port);

	buffer = buffer_create(1024, buffer_fixed, 1);
	
	auth(socket, buffer);
}