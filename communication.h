#include <sys/types.h>
#include <sys/socket.h>

#ifdef USE_UNIX_SOCKET
#include <sys/un.h>
#define SOCKET_TYPE AF_UNIX
#define SOCKET_ADDR struct sockaddr_un
#define SOCKET_ADDR_LEN sizeof(struct sockaddr_un)
#else
#include <netinet/in.h>
#include <arpa/inet.h>
#define SOCKET_TYPE AF_INET
#define SOCKET_ADDR struct sockaddr_in
#define SOCKET_ADDR_LEN sizeof(struct sockaddr_in)
#define PORT 6969
#endif

#define SOCKET_NAME "/tmp/pimc.sock"
#define HANDSHAKE_MSG "START"

typedef enum {
    SOCKOP_SUCCESS = 0,
    SOCKOP_ERROR,
    SOCKOP_DISCONNECTED,
} SocketOpResult;

SocketOpResult send_trace(int sockfd, EnergyTrace tr);
SocketOpResult recv_trace(int sockfd, EnergyTrace *tr);

int start_client();
int start_server();

SocketOpResult send_chars(int sockfd, const char *msg);
SocketOpResult recv_chars(int sockfd, char *buffer);

#ifdef COMMUNICATION_IMPLEMENTATION

SocketOpResult send_trace(int sockfd, EnergyTrace tr)
{
    if (send(sockfd, &tr.count, sizeof(tr.count), 0) != sizeof(tr.count)) {
        perror("send");
        return SOCKOP_ERROR; 
    }

    ssize_t sent = send(sockfd, tr.items, tr.count*sizeof(tr.items[0]), 0);
    if (sent != (ssize_t) (tr.count*sizeof(tr.items[0]) )) {
        perror("send");
        return SOCKOP_ERROR; 
    } 

    return SOCKOP_SUCCESS; 
}

SocketOpResult recv_trace(int sockfd, EnergyTrace *tr)
{
    memset(tr, 0, sizeof(EnergyTrace)); 

    ssize_t ret = recv(sockfd, &tr->count, sizeof(tr->count), 0);
    if (ret == 0) {
        return SOCKOP_DISCONNECTED; 
    }
    
    tr->items = arena_alloc(&arena, tr->count*sizeof(double));
    memset(tr->items, 0, tr->count*sizeof(double));
    tr->capacity = tr->count;

    ssize_t r = recv(sockfd, tr->items, sizeof(double)*tr->count, 0);
    if (r != (ssize_t) (tr->count*sizeof(double))) {
        perror("recv");
        return SOCKOP_ERROR; 
    } 

    return SOCKOP_SUCCESS; 
}

int start_client()
{
    int sockfd = socket(SOCKET_TYPE, SOCK_STREAM, 0);

#ifdef USE_UNIX_SOCKET
    SOCKET_ADDR server_addr;
    server_addr.sun_family = SOCKET_TYPE; 
    strcpy(server_addr.sun_path, SOCKET_NAME);
#else
    SOCKET_ADDR server_addr;
    server_addr.sin_family = SOCKET_TYPE;
    server_addr.sin_port = htons(PORT);

    // set ip in binary form
    if (inet_pton(AF_INET, "127.0.0.1", &server_addr.sin_addr) <= 0) {
        perror("inet_pton");
        exit(1);
    }
#endif

    if (connect(sockfd, (struct sockaddr *) &server_addr, sizeof(server_addr)) < 0) {
        perror("connect");
        return -1; 
    }

    if (send_chars(sockfd, HANDSHAKE_MSG) < 0) return -1;
        
    char *msg = NULL;
    SocketOpResult r = recv_chars(sockfd, msg);
    if (r != SOCKOP_SUCCESS) {
        return -1;
    }

    if (msg && (strcmp(msg, HANDSHAKE_MSG) != 0)) {
        fprintf(stderr, "[client] Connection is not established\n");
        return -1; 
    }

    return sockfd;
}

int start_server()
// TODO: maybe use non-blocking socket connection to postpone handling the message
{
    int server_socket = socket(SOCKET_TYPE, SOCK_STREAM, 0);
    if (server_socket < 0) {
        fprintf(stderr, "ERROR: could not create the server socket\n");
        exit(1);
    }

#ifdef USE_UNIX_SOCKET 
    unlink(SOCKET_NAME);
    
    SOCKET_ADDR server_addr;
    memset(&server_addr, 0, SOCKET_ADDR_LEN);
    server_addr.sun_family = SOCKET_TYPE; 
    strcpy(server_addr.sun_path, SOCKET_NAME);
#else
    SOCKET_ADDR server_addr;
    memset(&server_addr, 0, SOCKET_ADDR_LEN);
    server_addr.sin_family = SOCKET_TYPE;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    server_addr.sin_port = htons(PORT);
#endif 


    int ret;
    ret = bind(server_socket, (struct sockaddr *) &server_addr, sizeof(server_addr));
    if (ret < 0) {
        perror("bind");
        exit(1);
    }

    // 5 - connection queues?
    ret = listen(server_socket, 5);     
    if (ret < 0) {
        perror("listen");
        exit(1);
    }
    
    int data_socket;
 
    SOCKET_ADDR client_addr;
    int clen = SOCKET_ADDR_LEN; 

    data_socket = accept(server_socket, (struct sockaddr *) &client_addr, (socklen_t*) &clen);
    if (data_socket < 0) {
        perror("accept");
        exit(1);
    }

#ifndef USE_UNIX_SOCKET
    printf("Connection accepted from client IP address %s and %d\n",
            inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
#endif

    char *msg = NULL;
    do {
        SocketOpResult r = recv_chars(data_socket, msg);
        switch (r) {
            case SOCKOP_ERROR: exit(1);
            case SOCKOP_DISCONNECTED: exit(1); 
            case SOCKOP_SUCCESS:
        }
    } while (msg && (strcmp(msg, HANDSHAKE_MSG) != 0));
    printf("server: received hanshake\n");

    if (send_chars(data_socket, HANDSHAKE_MSG) < 0) return -1; 

    return data_socket;
}

SocketOpResult send_chars(int sockfd, const char *msg)
// prepend the char buffer with a "header" 
{
    uint32_t msg_length = strlen(msg);
    if (send(sockfd, &msg_length, sizeof(msg_length), 0) != sizeof(msg_length)) {
        perror("send");
        return SOCKOP_ERROR; 
    }
    
    ssize_t sent = send(sockfd, msg, msg_length, 0);
    if (sent != (ssize_t) msg_length) {
        perror("send");
        return SOCKOP_ERROR; 
    }
    
    return SOCKOP_SUCCESS;
}

SocketOpResult recv_chars(int sockfd, char *buffer)
{
    // assume that the message is prepended with a "header" that contains the length of the message 
    uint32_t msg_length = 0;
    ssize_t ret = recv(sockfd, &msg_length, sizeof(msg_length), 0);
    if (ret == 0) {
        return SOCKOP_DISCONNECTED; 
    }

    buffer = arena_alloc(&arena, msg_length*sizeof(char));
    memset(buffer, 0, msg_length*sizeof(char));

    ssize_t r = recv(sockfd, buffer, sizeof(char)*msg_length, 0);
    if (r != msg_length) {
        perror("recv");
        return SOCKOP_ERROR;
    } 

    return SOCKOP_SUCCESS; 
} 
#endif // COMMUNICATION_IMPLEMENTATION

