#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <errno.h>
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
} ConnectionStatus;

typedef enum {
    MSG_CHAR = 0,
    MSG_INT,
    MSG_DOUBLE, 
} MessageKind;

typedef struct {
    MessageKind kind;
    int length;
} MessageHeader;  

// Keep in mind 'flexible array member'?
typedef struct {
    uint32_t bytes_length; // 4 bytes
    MessageKind kind;      // 4 bytes
    void *payload;         // 8 bytes
} Message;


// TODO: send named messages: name+value
//       so that we could assert that the intended value has been received
//       On a server side we could put them into a map 
//
// Also we could built up the 'state' (#clients, beta etc.) and send it as a single structure,
// but this would mean that we should change the code for a server if a different client would want to add another field in the 'state' structure.
// So named parameter approach seems to be better.. 

#define HEADER_SIZE sizeof(MessageHeader)

static bool _verbose = true;
static char *_client_ip = NULL; // deallocation?

static int server_socket = 0;

int initClient();

void initServer();
ConnectionStatus acceptClientConnection(int *data_socket);

void set_verbose(bool flag); 
char *get_client_ip();

// serialize: message -> bytes stream
// deserialize: bytes stream -> message
void serialize_message(Message *message, uint8_t **bytes, uint32_t *bytes_length);
Message deserialize_message(uint8_t *bytes, uint32_t bytes_length);

// TODO: save sockfd in the local state?
SocketOpResult sendChars(int sockfd, const char *msg);
SocketOpResult recvChars(int sockfd, char *buffer);

SocketOpResult sendDoubleArray(int sockfd, double *data, size_t count);
SocketOpResult recvDoubleArray(int sockfd, double **data, size_t *count);

SocketOpResult sendInt(int sockfd, int value);
SocketOpResult recvInt(int sockfd, int *value);

SocketOpResult sendDouble(int sockfd, double value);
SocketOpResult recvDouble(int sockfd, double *value);


#ifdef PROTOCOL_IMPLEMENTATION 

void set_verbose(bool flag) {
    _verbose = flag;
}

char *get_client_ip() {
    assert(_client_ip != NULL);
    return _client_ip;
}

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

void initServer()
{
    server_socket = socket(SOCKET_TYPE, SOCK_STREAM, 0);
    if (server_socket < 0) {
        fprintf(stderr, "ERROR: could not create the server socket\n");
        exit(1);
    }

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
    do {
        SocketOpResult r = recvChars(*data_socket, msg);
        switch (r) {
            case SOCKOP_ERROR: exit(1);
            case SOCKOP_DISCONNECTED: exit(1); 
            case SOCKOP_SUCCESS: break;
            case SOCKOP_WAITING: assert(false);
        }
    } while (msg && (strcmp(msg, HANDSHAKE_MSG) != 0));
    printf("server: received hanshake\n");

    if (sendChars(*data_socket, HANDSHAKE_MSG) < 0) return -1; 

    return CONNECTION_ESTABLISHED; 
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
    if (bytes_recv == -1) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
            return SOCKOP_WAITING; 
        } else {
            perror("recv");
            return SOCKOP_ERROR;
        }
    } 

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

