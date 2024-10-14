#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <libssh2.h>
#include <lissh2_publickey.h>
//#include <libssh2_sftp.h>

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

#define PROTOCOL_IMPLEMENTATION
#include "protocol.h"

#define BUFFER_SIZE 8192 

int main()
{
    // local end of the tunnel 
    int local_socket; 
    if ((local_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
        perror("socket failed");
        exit(1);
    }

    int remote_socket;
    if ((remote_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
        perror("socket failed");
        exit(1);
    }

    struct sockaddr_in address;
    int addrlen = sizeof(address);
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(22);
    
    inet_pton(AF_INET, "10.62.5.57", &address.sin_addr);
    if (connect(remote_socket, (struct sockaddr*) address, addrlen) < 0) {
        perror("connect");
        exit(1);
    }

    // Initialize the SSH session
    libssh2_session *session = libssh2_session_init();
    if (session == NULL) {
        perror("libssh2_session_init");
        exit(1);
    }

    // Authenticate with the SSH server
    if (libssh2_userauth_password(session, "a.finenko", "1hpTPyvJew6v") != 0) {
        perror("libssh2_userauth_password");
        exit(1);
    }

    // create a channel for the SSH connection
    libssh2_channel *channel = libssh2_channel_open_session(session);
    if (channel == NULL) {
        perror("libssh2_channel_open_session failed");
        exit(1);
    }

    // Create a channel for the forwarded connection
    libssh2_channel *forward_channel = libssh2_channel_forward_listen_ex(session, "localhost", PORT, &local_socket);
    if (forward_channel == NULL) {
        perror("libssh2_channel_forward_listen_ex failed");
        exit(1);
    }

    if (bind(server_fd, (struct sockaddr*) &address, addrlen) < 0) {
        perror("bind failed");
        exit(1);
    }

    if (listen(server_fd, 3) < 0) {
        perror("listen failed");
        exit(1);
    }

    printf("Server listening on port %d...\n", PORT);

    int data_socket;
    if ((data_socket = accept(server_fd, (struct sockaddr*) &address, (socklen_t*) &addrlen)) < 0) {
        perror("accept failed");
        exit(1);
    }

    printf("Connection accepted\n");

    char response[BUFFER_SIZE] = {0};
    void* buffer = malloc(BUFFER_SIZE*sizeof(char));
    assert(buffer != NULL);

    ssize_t bytes_recv;

    while (1) {
        MessageHeader header; 
        bytes_recv = recv(data_socket, &header, HEADER_SIZE, 0);
        assert(bytes_recv == HEADER_SIZE);
        
        memset(buffer, 0, BUFFER_SIZE*sizeof(char));

        if (header.kind == MSG_CHAR) {
            bytes_recv = recv(data_socket, buffer, header.length*sizeof(char), 0);
            if (bytes_recv != (ssize_t) (header.length*sizeof(char))) {
                perror("recv");
                exit(1); 
            } 
            printf("Received text message: %s\n", (char*) buffer);

        } else if (header.kind == MSG_DOUBLE) {
            bytes_recv = recv(data_socket, buffer, header.length*sizeof(double), 0);
            if (bytes_recv != (ssize_t) (header.length*sizeof(double))) {
                exit(1);
            }
            printf("Received %d doubles: ", header.length);

            for (int i = 0; i < header.length; ++i) {
                printf("%.3f ", ((double*) buffer)[i]);
            }
            printf("\n");
        } else {
            bytes_recv = recv(data_socket, buffer, header.length*sizeof(int), 0);
            if (bytes_recv != (ssize_t) (header.length*sizeof(int))) {
                exit(1);
            }
            printf("Received %d integers: ", header.length);

            for (int i = 0; i < header.length; ++i) {
                printf("%d ", ((int*) buffer)[i]);
            }
            printf("\n");
        }

        printf("Enter response: ");
        char *ret = fgets(response, BUFFER_SIZE, stdin);
        assert(ret != NULL);

        response[strcspn(response, "\n")] = 0;

        if (strlen(response) > 0) {
            sendChars(data_socket, response);
            printf("response sent\n");
        } else {
            printf("skipping the response\n");
        }
    }

    return 0;
}
