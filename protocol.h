#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <errno.h>

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
#endif

#define SOCKET_NAME "/tmp/pimc.sock"
#define HANDSHAKE_MSG "START"

#define SERVER_IP "127.0.0.1"
#define PORT 6969

typedef enum {
    SOCKOP_SUCCESS = 0,
    SOCKOP_ERROR,
    SOCKOP_DISCONNECTED,
} SocketOpResult;

typedef enum {
    MSG_CHAR = 0,
    MSG_INT,
    MSG_DOUBLE, 
} MessageKind;

typedef struct {
    MessageKind kind;
    int length;
} MessageHeader;  

#define HEADER_SIZE sizeof(MessageHeader)

static bool _verbose = true;

// TODO: pass an arena through init functions? 
int initClient();
int initServer(bool verbose);

// TODO: save sockfd in the local state?
SocketOpResult sendChars(int sockfd, const char *msg);
SocketOpResult recvChars(int sockfd, char *buffer);

SocketOpResult sendDoubleArray(int sockfd, double *data, size_t count);
SocketOpResult recvDoubleArray(int sockfd, double **data, size_t *count);

SocketOpResult sendInt(int sockfd, int value);
SocketOpResult recvInt(int sockfd, int *value);

SocketOpResult sendDouble(int sockfd, double value);
SocketOpResult recvDouble(int sockfd, double *value);

// TODO: send named parameter? string+value

#ifdef PROTOCOL_IMPLEMENTATION 

/*
SocketOpResult recv_trace(int sockfd, EnergyTrace *tr)
{
    memset(tr, 0, sizeof(EnergyTrace)); 

    size_t count = 0;
    ssize_t ret = recv(sockfd, &count, sizeof(count), 0);
    if (ret == 0) {
        return SOCKOP_DISCONNECTED;
    }
    if (errno > 0) {
        printf("Error on receving trace size: %s\n", strerror(errno));
        return SOCKOP_ERROR;
    } 
    
    printf("Requesting to alloc %zu bytes\n", count);
    tr->items = arena_alloc(&arena, count*sizeof(double));
    memset(tr->items, 0, count*sizeof(double));

    ssize_t r = recv(sockfd, tr->items, sizeof(double)*count, 0);
    if (r != (ssize_t) (count*sizeof(double))) {
        perror("recv");
        return SOCKOP_ERROR; 
    }

    if (errno > 0) {
        printf("Error on receiving trace: %s\n", strerror(errno));
    }

    tr->count = count; 
    tr->capacity = count;
    printf("Successfully received %zu elements of trace\n", tr->count);

    return SOCKOP_SUCCESS; 
}
*/

int initClient()
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
    if (inet_pton(AF_INET, SERVER_IP, &server_addr.sin_addr) <= 0) {
        perror("inet_pton");
        exit(1);
    }
#endif

    if (connect(sockfd, (struct sockaddr *) &server_addr, sizeof(server_addr)) < 0) {
        perror("connect");
        return -1; 
    }

    if (sendChars(sockfd, HANDSHAKE_MSG) < 0) return -1;
        
    char *msg = NULL;
    SocketOpResult r = recvChars(sockfd, msg);
    if (r != SOCKOP_SUCCESS) {
        return -1;
    }

    if (msg && (strcmp(msg, HANDSHAKE_MSG) != 0)) {
        fprintf(stderr, "[client] Connection is not established\n");
        return -1; 
    }

    return sockfd;
}

int initServer(bool verbose)
// TODO: maybe use non-blocking socket connection to postpone handling the message
{
    _verbose = verbose;

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
   
    if (_verbose) { 
        printf("Server listening on port %d...\n", PORT);
    }
 
    SOCKET_ADDR client_addr;
    int clen = SOCKET_ADDR_LEN; 
    int data_socket;

    data_socket = accept(server_socket, (struct sockaddr *) &client_addr, (socklen_t*) &clen);
    if (data_socket < 0) {
        perror("accept");
        exit(1);
    }

#ifndef USE_UNIX_SOCKET
    if (_verbose) {
        printf("Connection accepted from client IP address %s and %d\n",
                inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
    }
#endif

    char *msg = NULL;
    do {
        SocketOpResult r = recvChars(data_socket, msg);
        switch (r) {
            case SOCKOP_ERROR: exit(1);
            case SOCKOP_DISCONNECTED: exit(1); 
            case SOCKOP_SUCCESS: break;
        }
    } while (msg && (strcmp(msg, HANDSHAKE_MSG) != 0));
    printf("server: received hanshake\n");

    if (sendChars(data_socket, HANDSHAKE_MSG) < 0) return -1; 

    return data_socket;
}

SocketOpResult sendChars(int sockfd, const char *msg)
{
    ssize_t bytes_sent;
    uint32_t msg_length = strlen(msg);

    MessageHeader header = {
        .kind = MSG_CHAR,
        .length = msg_length, 
    };

    bytes_sent = send(sockfd, &header, 1*HEADER_SIZE, 0);
    if (bytes_sent != (ssize_t) HEADER_SIZE) {
        perror("send");
        return SOCKOP_ERROR;
    } 

    bytes_sent = send(sockfd, msg, msg_length*sizeof(char), 0);
    if (bytes_sent != (ssize_t) msg_length) {
        perror("send");
        return SOCKOP_ERROR; 
    }
    
    return SOCKOP_SUCCESS;
}

SocketOpResult recvChars(int sockfd, char *buffer)
{
    MessageHeader header = {0};
    
    ssize_t bytes_recv = recv(sockfd, &header, 1*HEADER_SIZE, 0);
    assert(bytes_recv == HEADER_SIZE);
    assert(header.kind == MSG_CHAR); 

    buffer = arena_alloc(&arena, header.length*sizeof(char));
    memset(buffer, 0, header.length*sizeof(char));

    bytes_recv = recv(sockfd, buffer, header.length*sizeof(char), 0); 
    if (bytes_recv != (ssize_t) (header.length*sizeof(char))) {
        perror("recv");
        return SOCKOP_ERROR;
    } 

    return SOCKOP_SUCCESS; 
} 

// TODO: refactor sending the header
SocketOpResult sendDouble(int sockfd, double value)
{
    ssize_t bytes_sent;
    
    MessageHeader header = {
        .kind = MSG_DOUBLE,
        .length = 1,
    };

    bytes_sent = send(sockfd, &header, 1*HEADER_SIZE, 0);
    if (bytes_sent != (ssize_t) HEADER_SIZE) {
        perror("send");
        return SOCKOP_ERROR;
    }

    bytes_sent = send(sockfd, &value, 1*sizeof(double), 0);
    if (bytes_sent != sizeof(double)) {
        perror("send");
        return SOCKOP_ERROR;
    } 

    return SOCKOP_SUCCESS;
}

SocketOpResult recvDouble(int sockfd, double *buffer)
{
    MessageHeader header = {0};
    
    ssize_t bytes_recv = recv(sockfd, &header, 1*HEADER_SIZE, 0);
    assert(bytes_recv == HEADER_SIZE);
    assert(header.kind == MSG_DOUBLE); 

    bytes_recv = recv(sockfd, buffer, 1*sizeof(double), 0); 
    if (bytes_recv != 1*sizeof(double)) {
        perror("recv");
        return SOCKOP_ERROR;
    }

    if (_verbose) {
        printf("success: received double %.3lf\n", *buffer);
    } 

    return SOCKOP_SUCCESS; 
} 

SocketOpResult sendInt(int sockfd, int value)
{
    ssize_t bytes_sent;
    
    MessageHeader header = {
        .kind = MSG_INT,
        .length = 1,
    };

    bytes_sent = send(sockfd, &header, 1*HEADER_SIZE, 0);
    if (bytes_sent != (ssize_t) HEADER_SIZE) {
        perror("send");
        return SOCKOP_ERROR;
    }

    bytes_sent = send(sockfd, &value, 1*sizeof(int), 0);
    if (bytes_sent != sizeof(int)) {
        perror("send");
        return SOCKOP_ERROR;
    } 

    return SOCKOP_SUCCESS;
}

SocketOpResult recvInt(int sockfd, int *buffer)
{
    MessageHeader header = {0};
    
    ssize_t bytes_recv = recv(sockfd, &header, 1*HEADER_SIZE, 0);
    assert(bytes_recv == HEADER_SIZE);
    assert(header.kind == MSG_INT); 

    bytes_recv = recv(sockfd, buffer, 1*sizeof(int), 0); 
    if (bytes_recv != 1*sizeof(int)) {
        perror("recv");
        return SOCKOP_ERROR;
    } 
    
    if (_verbose) {
        printf("success: received int %d\n", *buffer);
    } 

    return SOCKOP_SUCCESS; 
} 

SocketOpResult sendDoubleArray(int sockfd, double *data, size_t count)
{
    ssize_t bytes_sent;
    
    MessageHeader header = {
        .kind = MSG_DOUBLE,
        .length = (int) count,
    };

    bytes_sent = send(sockfd, &header, 1*HEADER_SIZE, 0);
    if (bytes_sent != (ssize_t) HEADER_SIZE) {
        perror("send");
        return SOCKOP_ERROR;
    }

    bytes_sent = send(sockfd, data, count*sizeof(double), 0);
    if (bytes_sent != (ssize_t) (count*sizeof(double))) {
        perror("send");
        return SOCKOP_ERROR;
    } 

    return SOCKOP_SUCCESS;
}

SocketOpResult recvDoubleArray(int sockfd, double **data, size_t *count)
{
    MessageHeader header = {0};
    
    ssize_t bytes_recv = recv(sockfd, &header, 1*HEADER_SIZE, 0);
    assert(bytes_recv == HEADER_SIZE);
    assert(header.kind == MSG_DOUBLE); 

    *data = arena_alloc(&arena, header.length*sizeof(double));
    memset(*data, 0, header.length*sizeof(double));

    bytes_recv = recv(sockfd, *data, header.length*sizeof(double), 0); 
    if (bytes_recv != (ssize_t) (header.length*sizeof(double))) {
        perror("recv");
        return SOCKOP_ERROR;
    } 

    *count = header.length;

    if (_verbose) {
        printf("success: received array of double of %d elements: %.3lf ...\n", header.length, (*data)[0]);
    }

    return SOCKOP_SUCCESS; 
} 

#endif // PROTOCOL_IMPLEMENTATION 

