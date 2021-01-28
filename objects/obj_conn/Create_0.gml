address = "192.168.43.78";
port = 6621;

socket = network_create_socket(network_socket_tcp);
network_connect(socket, address, port);

buffer = buffer_create(1024, buffer_fixed, 0);