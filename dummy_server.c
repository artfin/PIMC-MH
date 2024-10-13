#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

#define PROTOCOL_IMPLEMENTATION
#include "protocol.h"

#define BUFFER_SIZE 8192 

int main()
{
    int server_fd;
    struct sockaddr_in address;
    int addrlen = sizeof(address);
    
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
        perror("socket failed");
        exit(1);
    }

    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT);

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
