#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>

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
#define MAX_CONNECTIONS 5

#define MAX_NAME_SIZE 20

typedef enum {
    SOCKOP_SUCCESS = 0,
    SOCKOP_WAITING, 
    SOCKOP_ERROR,
    SOCKOP_DISCONNECTED,
} SocketOpResult;

typedef enum {
    NO_CONNECTION,
    CONNECTION_ESTABLISHED, 
    CONNECTION_ERROR, // TODO: we may need to return from the 'accept' function to clean up the broken connection
    DISCONNECTED,
} ConnectionStatus;

typedef enum {
    MSG_CHAR = 0,
    MSG_INT32,
    MSG_FLOAT64,
    MSG_NAMED_FLOAT64,
} MessageKind;


// TODO: send named messages: name+value
//       so that we could assert that the intended value has been received
//       On a server side we could put them into a map 
//
// Also we could built up the 'state' (#clients, beta etc.) and send it as a single structure,
// but this would mean that we should change the code for a server if a different client would want to add another field in the 'state' structure.
// So named parameter approach seems to be better.. 

#define assert_sockop_result(r)                                             \
    do {                                                                    \
        if (r != SOCKOP_SUCCESS) {                                          \
            printf("%s:%d socket operation failed\n", __FILE__, __LINE__);  \
            return_defer(1);                                                \
        }                                                                   \
    } while(0)                                                              \

static bool _verbose = true;
static char *_client_ip = NULL; // deallocation?

// --------------------
static int server_socket = 0;

ConnectionStatus acceptClientConnection(int *data_socket);
void initServer();

char *get_client_ip();
// --------------------

// --------------------
int initClient();
// --------------------

void set_verbose(bool flag); 

// sizeof(Message) = 8
typedef struct {
    uint32_t size;      // [4 bytes] the length in bytes of [size + kind + payload] = 8 + [payload] 
    MessageKind kind;   // [4 bytes]
    uint8_t payload[];  // [0 bytes (if empty)]
} __attribute__((packed)) Message;

// serialize: message -> bytes stream
// deserialize: bytes stream -> message
// void serialize_message(Message *message, uint32_t *bytes_length, uint8_t **bytes); // TODO: what is the approach here? 
Message* deserialize_message(MessageKind kind, uint32_t bytes_length, uint8_t *bytes);

SocketOpResult sendNamedFloat64Array(int sockfd, const char *name, double *data, size_t count);
SocketOpResult recvNamedFloat64Array(int sockfd, char **name, double **data, size_t *count);

SocketOpResult sendFloat64Array(int sockfd, double *data, size_t count);
SocketOpResult recvFloat64Array(int sockfd, double **data, size_t *count);

SocketOpResult sendInt32(int sockfd, int32_t value);
SocketOpResult recvInt32(int sockfd, int32_t *value);

SocketOpResult sendFloat64(int sockfd, double value);
SocketOpResult recvFloat64(int sockfd, double *value);

SocketOpResult sendFixedLengthString(int sockfd, const char *buffer);
SocketOpResult recvFixedLengthString(int sockfd, char **buffer, size_t *count);

#ifdef PROTOCOL_IMPLEMENTATION 

void set_verbose(bool flag) {
    _verbose = flag;
}

char *get_client_ip() {
    assert(_client_ip != NULL);
    return _client_ip;
}

Message* deserialize_message(MessageKind kind, uint32_t payload_length, uint8_t *payload)
{
    Message *message = (Message*) arena_alloc(&arena, sizeof(message) + payload_length); 
    message->size = sizeof(message) + payload_length;
    message->kind = kind;
    memcpy(message->payload, payload, payload_length); 

    return message; 
}


SocketOpResult w_send(int sockfd, const void *buf, size_t len) {
    ssize_t bytes_sent = send(sockfd, buf, len, 0);
    if (bytes_sent == 0) {
        return SOCKOP_DISCONNECTED;
    } else if (bytes_sent != (ssize_t) len) {
        perror("send");
        return SOCKOP_ERROR;
    }

    return SOCKOP_SUCCESS;
} 

SocketOpResult w_recv(int sockfd, void *buf, size_t len) {
    ssize_t bytes_recv = recv(sockfd, buf, len, 0);
    if (bytes_recv == 0) {
        return SOCKOP_DISCONNECTED;
    } else if (bytes_recv == -1) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
            return SOCKOP_WAITING; 
        } else {
            perror("recv");
            return SOCKOP_ERROR;
        }
    }
    if (bytes_recv != (ssize_t) len) {
        return SOCKOP_ERROR;
    }

    return SOCKOP_SUCCESS;
}

SocketOpResult sendInt32(int sockfd, int32_t value)
{
    Message *msg = deserialize_message(MSG_INT32, sizeof(int), (uint8_t*) &value);
    
    if (_verbose) {
        printf("(sendInt32): message contents = [int] %d\n", value);
    } 
    
    return w_send(sockfd, msg, msg->size);
}

SocketOpResult recvInt32(int sockfd, int32_t *value)
{
    size_t sz = sizeof(Message) + sizeof(int);
    
    void *bytes = arena_alloc(&arena, sz);
    memset(bytes, 0, sz);
    
    SocketOpResult res = w_recv(sockfd, bytes, sz);
    if (res != SOCKOP_SUCCESS) return res; 

    Message *msg = (Message*) bytes;
    assert(msg->kind == MSG_INT32);
    memcpy(value, msg->payload, sizeof(int));
    
    if (_verbose) {
        printf("(recvInt32): message contents = [int] %d\n", *value);
    } 

    return SOCKOP_SUCCESS; 
} 

SocketOpResult sendFloat64(int sockfd, double value)
{
    Message *msg = deserialize_message(MSG_FLOAT64, sizeof(double), (uint8_t*) &value); 
    return w_send(sockfd, msg, msg->size); 
}

SocketOpResult recvFloat64(int sockfd, double *value)
{
    size_t sz = sizeof(Message) + sizeof(double);

    void *bytes = arena_alloc(&arena, sz);
    memset(bytes, 0, sz);
    
    SocketOpResult res = w_recv(sockfd, bytes, sz);
    if (res != SOCKOP_SUCCESS) return res; 

    Message *msg = (Message*) bytes;
    assert(msg->kind == MSG_FLOAT64);
    memcpy(value, msg->payload, sizeof(double));

    if (_verbose) {
        printf("(recvFloat64): message contents = [double] %.3lf\n", *value);
    } 

    return SOCKOP_SUCCESS; 
} 

SocketOpResult sendFloat64Array(int sockfd, double *data, size_t count)
{
    Message *msg = deserialize_message(MSG_FLOAT64, count*sizeof(double), (uint8_t*) data);
    return w_send(sockfd, msg, msg->size); 
}

SocketOpResult recvFloat64Array(int sockfd, double **data, size_t *count)
{
    SocketOpResult res;
    uint32_t sz;
    
    res = w_recv(sockfd, &sz, sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res; 
    
    void *bytes = arena_alloc(&arena, sz);
    memset(bytes, 0, sz);
    memcpy(bytes, &sz, sizeof(uint32_t));

    res = w_recv(sockfd, bytes + sizeof(uint32_t), sz - sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res; 
    
    Message *msg = (Message*) bytes;
    assert(msg->kind == MSG_FLOAT64);

    size_t payload_len = sz - sizeof(uint32_t) - sizeof(MessageKind); 
    assert(payload_len % sizeof(double) == 0);
    
    *count = payload_len / sizeof(double);
    *data = bytes + sizeof(uint32_t) + sizeof(MessageKind);

    if (_verbose) {
        printf("(recvFloat64Array): message contents = [%u bytes, %zu float64s] %.3lf ...\n", sz, *count, (*data)[0]);
    }

    return SOCKOP_SUCCESS; 
} 

SocketOpResult sendFixedLengthString(int sockfd, const char *buffer)
{
    Message *msg = deserialize_message(MSG_CHAR, strlen(buffer)*sizeof(char), (uint8_t*) buffer);
    return w_send(sockfd, msg, msg->size); 
}

SocketOpResult recvFixedLengthString(int sockfd, char **buffer, size_t *count)
{
    SocketOpResult res;
    uint32_t sz;
    
    res = w_recv(sockfd, &sz, sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res;

    void *bytes = arena_alloc(&arena, sz);
    memset(bytes, 0, sz);
    memcpy(bytes, &sz, sizeof(uint32_t));

    res = w_recv(sockfd, bytes + sizeof(uint32_t), sz - sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res; 
    
    Message *msg = (Message*) bytes;
    assert(msg->kind == MSG_CHAR);

    *buffer = bytes + sizeof(uint32_t) + sizeof(MessageKind);
    *count = sz - sizeof(uint32_t) - sizeof(MessageKind);  
 
    if (_verbose) {
        printf("(recvFixedLengthString): message contents = [%zu chars]: %.*s\n", *count, (int) *count, *buffer);
    }

    return SOCKOP_SUCCESS; 
} 

int initClient()
{
    int sockfd = socket(SOCKET_TYPE, SOCK_STREAM, 0);
    
    // allows reusing the socket address when the previous run of the application
    // exited with an error and did not properly closed the socket
    int option = 1;
    setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &option, sizeof(option));

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

    if (sendFixedLengthString(sockfd, HANDSHAKE_MSG) < 0) return -1;
        
    char *msg = NULL;
    size_t msg_length = 0;
    SocketOpResult r = recvFixedLengthString(sockfd, &msg, &msg_length);
    if (r != SOCKOP_SUCCESS) {
        return -1;
    }

    if (msg && (strcmp(msg, HANDSHAKE_MSG) != 0)) {
        fprintf(stderr, "[client] Connection is not established\n");
        return -1; 
    }

    return sockfd;
}

void initServer()
{
    server_socket = socket(SOCKET_TYPE, SOCK_STREAM, 0);
    if (server_socket < 0) {
        fprintf(stderr, "ERROR: could not create the server socket: %s\n", strerror(errno));
        exit(1);
    }

    // allows reusing the socket address when the previous run of the application
    // exited with an error and did not properly closed the socket
    int option = 1;
    setsockopt(server_socket, SOL_SOCKET, SO_REUSEADDR, &option, sizeof(option));

    int flags = fcntl(server_socket, F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(server_socket, F_SETFL, flags); 

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

    ret = listen(server_socket, MAX_CONNECTIONS);     
    if (ret < 0) {
        perror("listen");
        exit(1);
    }
   
    if (_verbose) { 
        printf("Server listening on port %d...\n", PORT);
    }
}

// TODO: maybe use non-blocking socket connection to postpone handling the message
ConnectionStatus acceptClientConnection(int *data_socket)
{
    SOCKET_ADDR client_addr;
    int clen = SOCKET_ADDR_LEN; 
    int client_fd = 0;

    client_fd = accept(server_socket, (struct sockaddr *) &client_addr, (socklen_t*) &clen);
    if (client_fd < 0) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
            return NO_CONNECTION; 
        } else {
            perror("accept");
            exit(1);
            // return CLIENT_CONNECTION_ERROR;
        }
    }

    *data_socket = client_fd;

#ifndef USE_UNIX_SOCKET
    if (_verbose) {
        printf("Connection accepted from client IP address %s:%d\n",
                inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
    }
    
    _client_ip = malloc(strlen(inet_ntoa(client_addr.sin_addr)));
    strcpy(_client_ip, inet_ntoa(client_addr.sin_addr));
#endif

    char *msg = NULL;
    size_t msg_length = 0;

    // do {
    // } while (strcmp(msg, HANDSHAKE_MSG) != 0);
    
    SocketOpResult r = recvFixedLengthString(*data_socket, &msg, &msg_length);
    if (r != SOCKOP_SUCCESS) return CONNECTION_ERROR;
    assert(strncmp(msg, HANDSHAKE_MSG, strlen(HANDSHAKE_MSG)) == 0);

    printf("server: received hanshake\n");
    if (sendFixedLengthString(*data_socket, HANDSHAKE_MSG) < 0) return -1; 

    return CONNECTION_ESTABLISHED; 
}

SocketOpResult sendNamedFloat64Array(int sockfd, const char *name, double *data, size_t count)
{
    uint32_t payload_length = MAX_NAME_SIZE*sizeof(char) + count*sizeof(double);

    Message *message = (Message*) arena_alloc(&arena, sizeof(message) + payload_length);
    message->size = sizeof(message) + payload_length;
    message->kind = MSG_NAMED_FLOAT64;

    assert(strlen(name) < MAX_NAME_SIZE);
    memcpy(message->payload, name, MAX_NAME_SIZE*sizeof(char));
    memcpy(message->payload + MAX_NAME_SIZE*sizeof(char), data, count*sizeof(double));
        
    return w_send(sockfd, message, message->size); 
}

SocketOpResult recvNamedFloat64Array(int sockfd, char **name, double **data, size_t *count)
{
    SocketOpResult res;
    uint32_t sz;
    
    res = w_recv(sockfd, &sz, sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res; 
    
    void *bytes = arena_alloc(&arena, sz);
    memset(bytes, 0, sz);
    memcpy(bytes, &sz, sizeof(uint32_t));

    res = w_recv(sockfd, bytes + sizeof(uint32_t), sz - sizeof(uint32_t));
    if (res != SOCKOP_SUCCESS) return res; 
    
    Message *msg = (Message*) bytes;
    assert(msg->kind == MSG_NAMED_FLOAT64);

    *name = bytes + sizeof(uint32_t) + sizeof(MessageKind);

    size_t payload_len = sz - sizeof(uint32_t) - sizeof(MessageKind) - MAX_NAME_SIZE*sizeof(char); 
    assert(payload_len % sizeof(double) == 0);
    
    *count = payload_len / sizeof(double);
    *data = bytes + sizeof(uint32_t) + sizeof(MessageKind) + MAX_NAME_SIZE*sizeof(char);

    if (_verbose) {
        printf("(recvNamedFloat64Array): message contents = [name = %.*s, %u bytes, %zu float64s] %.3lf ...\n", MAX_NAME_SIZE, *name, sz, *count, (*data)[0]);
    }
    
    return SOCKOP_SUCCESS; 
} 

#endif // PROTOCOL_IMPLEMENTATION 

