// Based on the notes by Adrian Del Maestro 
// "Path Integral Monte Carlo and the Worm algorithm in the Spatial Continuum"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>


#include <gsl/gsl_histogram.h>

#include "raylib.h"

#define FONT_SIZE_LOAD 160 
Font font = {0};

#define MT_GENERATE_CODE_IN_HEADER 0
#include "mtwist.h"
/*
 * double mt_drand(void)
 *   Return a pseudorandom double in [0,1) with 32 bits of randomness
 *
 * uint32_t mt_lrand(void);
 *   Generate 32-bit random value 
 *
 */
double generate_normal(double sigma) 
/*
 * Generate normally distributed variable using Box-Muller method
 */
{
    double U = mt_drand();
    double V = mt_drand();
    return sigma * sqrt(-2 * log(U)) * cos(2.0 * M_PI * V);
}

#define MAX_FLOAT_TEXTLEN 100

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

#define COMMON_IMPLEMENTATION
#include "common.h"

void subcmd_visualize(const char *program_path, int argc, char **argv);
void subcmd_run(const char *program_path, int argc, char **argv);
void subcmd_client(const char *program_path, int argc, char **argv);

Subcmd subcmds[] = {
    DEFINE_SUBCMD(visualize, "Visualize the chain modifications"),
    DEFINE_SUBCMD(run, "Run the PIMC calculation of the harmonic oscillator"),
    DEFINE_SUBCMD(client, "Start the client"),
};
#define SUBCMDS_COUNT (sizeof(subcmds)/sizeof(subcmds[0]))

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} EnergyTrace;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} PositionTrace;

typedef struct {
    EnergyTrace energies;
    PositionTrace positions;
} PIMC_Trace;

PIMC_Trace trace = {0};

// ---------------------------------------------------------
double lam = 0.5; // hbar^2/2m k_B
// ---------------------------------------------------------


// ---------------------------------------------------------
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

#define SOCKET_NAME "/tmp/pimc.sock"
#define HANDSHAKE_MSG "START"

typedef enum {
    SOCKOP_SUCCESS = 0,
    SOCKOP_ERROR,
    SOCKOP_DISCONNECTED,
} SocketOpResult;

SocketOpResult send_trace(int sockfd, EnergyTrace tr);
SocketOpResult recv_trace(int sockfd, EnergyTrace *tr);
// ---------------------------------------------------------

Path path = {0};

double V(double x) { return 0.5*x*x; }

double SHOExact(double beta) {
    //return 0.5/tanh(0.5/T);
    return 0.5/tanh(0.5*beta);
}

double PotentialAction(Path path, size_t tslice) 
{
    assert(tslice < path.numTimeSlices);

    double pot = 0.0;
    for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
        pot = pot + V(path.beads[tslice][ptcl]);
    }

    return path.tau * pot;
}

double PotentialEnergy(Path path) {
    double result = 0.0;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
            double x = path.beads[tslice][ptcl];
            result = result + V(x);
        }
    }

    return result / path.numTimeSlices;
}

double KineticEnergy(Path path) 
{
    double tot = 0.0;
    double norm = 1.0/(4.0*lam * path.tau*path.tau);

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        size_t tslicep1 = (tslice + 1) % path.numTimeSlices;

        for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
            double dx = path.beads[tslicep1][ptcl] - path.beads[tslice][ptcl];
            tot = tot - norm*dx*dx;
        }
    }

    return 0.5*path.numParticles/path.tau + tot/path.numTimeSlices;
}

double Energy(Path path) {
    return PotentialEnergy(path) + KineticEnergy(path);
}

typedef struct {
    double mean;
    double std;
} Stats;

Stats getStats(EnergyTrace t) 
{
    Stats s = {0};

    for (size_t i = 0; i < t.count; ++i) {
        s.mean = s.mean + t.items[i];
    }

    s.mean /= t.count;

    double r = 0;
    for (size_t i = 0; i < t.count; ++i) {
        r = r + (t.items[i] - s.mean) * (t.items[i] - s.mean); 
    }

    r = r / (t.count - 1);
    s.std = sqrt(r);

    return s;
}

int compar(const void *a, const void *b) {
    double *x = (double *) a;
    double *y = (double *) b;
    if (*x < *y) { 
        return -1;
    } else if (*x > *y) {
        return 1; 
    }
    return 0;
}

Stats getStatsEx(EnergyTrace t, size_t binSize)
// TODO: make sure that binning actually improves the std 
{
    qsort(t.items, t.count, sizeof(t.items[0]), compar);    

    int numBins = ceil((float) t.count / binSize);
    printf("numBins: %d\n", numBins);

    double *bins = (double*) arena_alloc(&arena, numBins * sizeof(double));
    memset(bins, 0.0, numBins * sizeof(double)); 

    int c = 0; // bin cursor
    for (size_t i = 0; i < t.count; ++i) 
    {
        if ((i > 0) && (i % binSize == 0)) {
            bins[c] /= binSize;
            c++;
        }

        bins[c] = bins[c] + t.items[i];
    }

    if (t.count % binSize != 0) {
        bins[c] /= (t.count % binSize);
    }

    //printf("Bins: ");
    //for (int i = 0; i < numBins; ++i) {
    //    printf("%.2lf ", bins[i]);
    //}
    //printf("\n");

    Stats s = {0};

    for (int i = 0; i < numBins; ++i) {
        s.mean = s.mean + bins[i];
    }

    s.mean /= numBins;

    double r = 0;
    for (int i = 0; i < numBins; ++i) {
        r = r + (bins[i] - s.mean) * (bins[i] - s.mean); 
    }

    r = r / (numBins - 1);
    s.std = sqrt(r);

    return s; 
}


int COM_Move(Path path, size_t ptcl)
{
    assert(ptcl < path.numParticles);
    int result = 0;

    double delta = 0.75;
    double shift = delta*(-1.0 + 2.0*mt_drand());

    double oldAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldAction = oldAction + PotentialAction(path, tslice);
    } 
    
    double *oldbeads = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeads, 0.0, path.numTimeSlices * sizeof(double));

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldbeads[tslice] = path.beads[tslice][ptcl];
        path.beads[tslice][ptcl] = path.beads[tslice][ptcl] + shift;
    }

    double newAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        newAction = newAction + PotentialAction(path, tslice);
    }
    
    // accept the move, or reject and restore the bead positions
    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));

    if (u < alpha) {
        return_defer(1);
    } else {
        for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
            path.beads[tslice][ptcl] = oldbeads[tslice];
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;
}

int Staging_Move(Path path, size_t ptcl)
// http://link.aps.org/doi/10.1103/PhysRevB.31.4234
{
    size_t stage_len = path.numTimeSlices/2; 
    assert(stage_len < path.numTimeSlices);

    int result = 0;

    size_t alpha_start = mt_lrand() % path.numTimeSlices;
    size_t alpha_end = (alpha_start + stage_len) % path.numTimeSlices; 

    double *oldbeads = (double*) arena_alloc(&arena, (stage_len - 1)*sizeof(double));
    memset(oldbeads, 0.0, (stage_len - 1)*sizeof(double));
    
    double oldAction = 0.0;

    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        oldbeads[i - 1] = path.beads[tslice][ptcl];
        oldAction = oldAction + PotentialAction(path, tslice);
    }

    double newAction = 0.0;
    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        size_t tslicem1 = (tslice - 1) % path.numTimeSlices;

        double tau1 = (stage_len - i) * path.tau;
        double avex = (tau1*path.beads[tslicem1][ptcl] + path.tau*path.beads[alpha_end][ptcl])/(path.tau+tau1);
        double sigma2 = 2.0*lam / (1.0/path.tau + 1.0/tau1);

        path.beads[tslice][ptcl] = avex + sqrt(sigma2)*generate_normal(1.0);
        //printf("avex: %.5lf; sigma2: %.5lf => pos = %.5lf\n", avex, sigma2, path.beads[tslice][ptcl]);
        newAction = newAction + PotentialAction(path, tslice); 
    }

    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));
    
    if (u < alpha) {
        return_defer(1);
    } else {
        for (size_t i = 1; i < stage_len; ++i) {
            size_t tslice = (alpha_start + i) % path.numTimeSlices;
            path.beads[tslice][ptcl] = oldbeads[i - 1];
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;

}

void pimc_driver(Path path, size_t numSteps, int sockfd, bool collect_positions)
{
    AcceptanceRate acc = {0};

    size_t equilSkip = 20000;
    size_t observableSkip = 400;

    size_t send_size = 1000;

    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %zu\n", equilSkip);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

    for (size_t step = 0; step < numSteps; ++step) 
    {
        for (size_t i = 0; i < path.numParticles; ++i) {
            size_t ptcl = mt_lrand() % path.numParticles;
            acc.CenterOfMass += COM_Move(path, ptcl);
        }

        for (size_t i = 0; i < path.numParticles; ++i) {
            size_t ptcl = mt_lrand() % path.numParticles;
            acc.Staging += Staging_Move(path, ptcl);
        }

       if ((step % observableSkip == 0) && (step > equilSkip)) {
            da_append(&trace.energies, Energy(path));

            if (collect_positions) {
                assert(path.numParticles == 1);
                for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                    da_append(&trace.positions, path.beads[tslice][0]);
                }
            }

            if (sockfd <= 0) continue; 
            if (trace.energies.count == send_size) {
                // @TODO: what if send was unsuccessful?
                send_trace(sockfd, trace.energies); 
                trace.energies.count = 0;
            }
       }
    }

    printf("--------------------------------------\n");
    printf("Acceptance Ratios:\n");
    printf("Center of Mass : %.3lf\n", (double) acc.CenterOfMass/numSteps/path.numParticles);
    printf("Staging        : %.3lf\n", (double) acc.Staging/numSteps/path.numParticles);
    printf("--------------------------------------\n");
}

void alloc_beads(Path *path)
{
    path->beads = (double**) malloc(path->numTimeSlices * sizeof(double*));
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        path->beads[i] = (double*) malloc(path->numParticles * sizeof(double));
        memset(path->beads[i], 0.0, path->numParticles*sizeof(double));
    }
}

void dealloc_beads(Path *path)
{
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        free(path->beads[i]);
    }

    free(path->beads);
}

#define COLOR_BACKGROUND GetColor(0x181818FF)

void update_draw_frame()
{
    static int simulation_step = 0; 

    if (IsKeyPressed(KEY_SPACE)) {
        printf("Simulation step: %d\n", simulation_step);

        if (simulation_step % 2 == 0) {
            for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
                COM_Move(path, ptcl);
            }
        } else {
            for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
                Staging_Move(path, ptcl);
            }
        }

        simulation_step++;
    }

    BeginDrawing();

    ClearBackground(COLOR_BACKGROUND);
   
    Color colors[] = {
        GetColor(0xF2AF29FF),
        GetColor(0xB52A2AFF),
        GetColor(0x7DD181FF),
    };
    #define COLORS_COUNT (sizeof(colors)/sizeof(colors[0]))

    int screen_width = GetScreenWidth();
    int screen_height = GetScreenHeight();

    int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));

    Rectangle world = {
        .x = 0.1*screen_width, 
        .y = 0.1*screen_height, 
        .width = simul_sz, 
        .height = simul_sz
    };

    DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);
    double GRID_MIN = -2.5;
    double GRID_MAX =  2.5;   

    double tau = world.height / (path.numTimeSlices - 1);

    double prevx, prevy;
    for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
        assert(ptcl < COLORS_COUNT);
        Color c = colors[ptcl];

        for (size_t i = 0; i < path.numTimeSlices; ++i) {
            double bead = path.beads[i][ptcl];
            if ((bead < GRID_MIN) || (bead > GRID_MAX)) {
                printf("%.5lf\n", bead);
            }

            double screenx = world.x + (bead - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
            double screeny = world.y + tau * i; 
            DrawCircle(screenx, screeny, 6.0, c);

            if (i > 0) {
                DrawLine(screenx, screeny, prevx, prevy, c);
            }

            prevx = screenx;
            prevy = screeny;
        }
    }

    EndDrawing(); 
}

void subcmd_visualize(const char *program_path, int argc, char **argv)
{
    (void) program_path;
    (void) argc;
    (void) argv;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "PIMC");
    SetExitKey(KEY_Q);

    double T = 1.0; // K
    path.numParticles = 1;
    path.numTimeSlices = 16;
    path.tau = 1.0/(T * path.numTimeSlices);
    alloc_beads(&path);
    
    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < path.numParticles; ++j) {
            path.beads[i][j] = 0.5 * (-1.0 + 2.0*mt_drand()); 
        }
    }

    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update_draw_frame();
    }

    CloseWindow();
}

void load_resources()
{
    int fileSize = 0;
    unsigned char* fileData = LoadFileData("resources/Alegreya-Regular.ttf", &fileSize);

    font.baseSize = FONT_SIZE_LOAD;
    font.glyphCount = 95;
    font.glyphs = LoadFontData(fileData, fileSize, FONT_SIZE_LOAD, 0, 95, FONT_SDF);
    Image atlas = GenImageFontAtlas(font.glyphs, &font.recs, 95, FONT_SIZE_LOAD, 4, 0);
    font.texture = LoadTextureFromImage(atlas);

    UnloadImage(atlas);
    UnloadFileData(fileData);
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
        //fprintf(stderr, "[server] connection closed\n");
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
        fprintf(stderr, "[server] connection closed\n");
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

int start_server()
// TODO: maybe use non-blocking socket connection to postpone handling the message
{
    int server_socket = socket(AF_UNIX, SOCK_STREAM, 0);
    if (server_socket < 0) {
        fprintf(stderr, "ERROR: could not create the server socket\n");
        exit(1);
    }
    
    unlink(SOCKET_NAME);
    
    // bind the socket
    struct sockaddr_un server_addr;
    memset(&server_addr, 0, sizeof(struct sockaddr_un));

    server_addr.sun_family = AF_UNIX;
    strcpy(server_addr.sun_path, SOCKET_NAME);

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
    struct sockaddr_un client_addr;
    int clen = sizeof(client_addr);

    data_socket = accept(server_socket, (struct sockaddr *) &client_addr, (socklen_t*) &clen);
    if (data_socket < 0) {
        perror("accept");
        exit(1);
    }

    return data_socket;
}

gsl_histogram* gsl_histogram_extend(gsl_histogram* h)
{
    size_t nbins = h->n;
    double add_bins = 1;
    
    double xmin = h->range[0];
    double xmax = h->range[nbins];
    double dx = h->range[1] - h->range[0];
    
    double new_xmin = xmin - add_bins*dx; 
    double new_xmax = xmax + add_bins*dx; 
    
    gsl_histogram *new_h = gsl_histogram_alloc(nbins + 2*add_bins);
    gsl_histogram_set_ranges_uniform(new_h, new_xmin, new_xmax);
        
    size_t nc = add_bins; // cursor over the new histogram
    for (size_t i = 0; i < nbins; ++i) {
        new_h->bin[nc++] = gsl_histogram_get(h, i); 
    }

    gsl_histogram_free(h);

    return new_h;
}

void subcmd_client(const char *program_path, int argc, char **argv)
{
    (void) program_path;
    (void) argc;
    (void) argv;
    
    int server_socket = socket(AF_UNIX, SOCK_STREAM, 0);

    struct sockaddr_un server_addr;
    server_addr.sun_family = AF_UNIX;
    strcpy(server_addr.sun_path, SOCKET_NAME);

    if (connect(server_socket, (struct sockaddr *) &server_addr, sizeof(server_addr)) < 0) {
        perror("connect");
        return; 
    }

    if (send_chars(server_socket, HANDSHAKE_MSG) < 0) return;
        
    char *msg = NULL;
    SocketOpResult r = recv_chars(server_socket, msg);
    switch (r) {
        case SOCKOP_ERROR: exit(1);
        case SOCKOP_DISCONNECTED: return; 
        case SOCKOP_SUCCESS:
    } 

    if (msg && (strcmp(msg, HANDSHAKE_MSG) != 0)) {
        fprintf(stderr, "[client] Connection is not established\n");
        return; 
    }

    double beta;
    size_t numTimeSlices;
    recv(server_socket, &beta, sizeof(double), 0);
    recv(server_socket, &numTimeSlices, sizeof(size_t), 0);
    
    double en_exact = SHOExact(beta);
    EnergyTrace tr = {0};

    size_t nbins = 20;
    gsl_histogram *h = gsl_histogram_alloc(nbins);
    gsl_histogram_set_ranges_uniform(h, -0.5, 0.5);
    size_t samples_count = 0;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "Client");
    SetExitKey(KEY_Q);
   
    load_resources();
    int font_size = 24;

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(COLOR_BACKGROUND);
        
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();
    
        double mean = gsl_histogram_mean(h);
        double sigma = gsl_histogram_sigma(h);

        nbins = h->n;
        double ymax = gsl_histogram_max_val(h);
        double xmin = h->range[0];
        double xmax = h->range[nbins];
         
        int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));

        Rectangle world = {
            .x = 0.1*screen_width, 
            .y = 0.1*screen_height, 
            .width = simul_sz, 
            .height = simul_sz
        };
            
        double rect_width = world.width / nbins; 

        for (size_t i = 0; i < nbins; ++i) {
            double height = gsl_histogram_get(h, i)/ymax * world.height; 
            double factor = 0.9;

            Rectangle r = {
                .x = world.x + i*rect_width,
                .y = world.y + world.height - factor*height, 
                .width = rect_width,
                .height = factor*height, 
            };

            DrawRectangleLinesEx(r, 2.0, YELLOW);
        }
        
        DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);

        if (samples_count > 0) {
            Vector2 mean_st = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y,
            };
            Vector2 mean_end = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y + world.height,
            };
            DrawLineEx(mean_st, mean_end, 2.0, RED);

            const char *buffer = TextFormat("%.3lf", mean);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + (mean - xmin)/(xmax - xmin)*world.width - 0.25*text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, RED);
        }
        
        {
            const char *buffer = TextFormat("%.3lf", xmin);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
        {
            const char *buffer = TextFormat("%.3lf", xmax);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + world.width - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
       
        // Statistics 
        if (samples_count > 0) {
            Vector2 stats_pos = { 
                .x = world.x + world.width + 50, 
                .y = world.y 
            };

            const char *buffer;
            buffer = TextFormat("beta: %.2f", beta);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);
            
            stats_pos.y += font_size; 
            buffer = TextFormat("time slices: %zu", numTimeSlices);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);

            stats_pos.y += font_size; 
            buffer = TextFormat("Samples: %zu", samples_count);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);
            
            stats_pos.y += font_size; 
            buffer = TextFormat("Number of bins: %zu", nbins);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);

            stats_pos.y += font_size; 
            buffer = TextFormat("Mean energy: %.5f", mean);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);
            
            stats_pos.y = stats_pos.y + font_size;    
            buffer = TextFormat("Exact energy: %.5f", en_exact);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);

            double exp_error = sigma/sqrt(samples_count);
            stats_pos.y += font_size;
            buffer = TextFormat("Error estimate: %.5f", exp_error);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);

            double actual_error = mean - en_exact; 
            stats_pos.y = stats_pos.y + font_size;    
            buffer = TextFormat("Actual error: %.5f", actual_error);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);

            double rel_error = fabs(actual_error) / en_exact;
            stats_pos.y = stats_pos.y + font_size;    
            buffer = TextFormat("Relative error: %.2f%%", rel_error*100.0);
            DrawTextEx(font, buffer, stats_pos, font_size, 0, WHITE);
        }

        EndDrawing();
        
        // @TODO: recv in non-blocking mode to still have client in responsive state  
        SocketOpResult r = recv_trace(server_socket, &tr);
        switch (r) {
            case SOCKOP_DISCONNECTED: continue; 
            case SOCKOP_ERROR: break;
            case SOCKOP_SUCCESS: 
        };
     
        for (size_t i = 0; i < tr.count; ++i) {
            while ((tr.items[i] < h->range[0] || (tr.items[i] > h->range[h->n]))) {
                h = gsl_histogram_extend(h);
                printf("extending histogram: %.5f -- %.5f\n", h->range[0], h->range[h->n]);
            } 
            gsl_histogram_increment(h, tr.items[i]);
        }
        samples_count += tr.count;

        arena_reset(&arena); 
    }

    close(server_socket);
    arena_free(&arena);
    gsl_histogram_free(h);
}


void subcmd_run(const char *program_path, int argc, char **argv)
{
    (void) program_path;
  
    bool opt_server = false; 
    bool opt_position_histogram = false;

    if (argc > 0) {
        char *opt = shift(&argc, &argv);

        if (strcmp(opt, "--server") == 0) {
            fprintf(stdout ,"-- Start the server\n");
            opt_server = true;
        } else if (strcmp(opt, "--position_histogram") == 0) {
            fprintf(stdout, "-- Collecting positions in the histogram\n");
            opt_position_histogram = true; 
        } else {
            // TODO: show available options
            fprintf(stderr, "ERROR: unknown option passed to `run` subcommand: %s\n", opt);
            exit(1); 
        }
        
        // TODO: accepting only one option for now
        assert(argc == 0);
    }
    
    uint32_t seed = mt_goodseed();
    mt_seed32(seed);
    
    path.numParticles = 1;
    path.numTimeSlices = 128;
    path.beta = 10.0;
    path.tau = path.beta/path.numTimeSlices;

    int sockfd = 0;
    if (opt_server) {
        sockfd = start_server();
       
        char *msg = NULL;
        do {
            SocketOpResult r = recv_chars(sockfd, msg);
            switch (r) {
                case SOCKOP_ERROR: exit(1);
                case SOCKOP_DISCONNECTED: { 
                    sockfd = 0;
                    break;
                }
                case SOCKOP_SUCCESS:
            }
        } while (msg && (strcmp(msg, HANDSHAKE_MSG) != 0));

        if (send_chars(sockfd, HANDSHAKE_MSG) < 0) return; 
    
        send(sockfd, &path.beta, sizeof(double), 0);
        send(sockfd, &path.numTimeSlices, sizeof(size_t), 0);
    }

    alloc_beads(&path);

    printf("Simulation parameters:\n");
    printf("Number of Particles   = %zu\n", path.numParticles);
    printf("Number of Time Slices = %zu\n", path.numTimeSlices);
    printf("beta                  = %.3lf\n", path.beta);
    printf("tau                   = %.3lf\n", path.tau);

    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < path.numParticles; ++j) {
            path.beads[i][j] = 0.5 * (-1.0 + 2.0*mt_drand()); 
        }
    }
    
    size_t MC_steps = 100 * 1000 * 1000;
    pimc_driver(path, MC_steps, sockfd, true);

    int binSize = 500;
    Stats s = getStatsEx(trace.energies, binSize);
    printf("Collected %zu values\n", trace.energies.count); 
    
    // NOTE: the error of the mean is actually the samples std divided by the sqrt(number-of-samples)  
    double en_err = s.std/sqrt(trace.energies.count); 
    printf("(PIMC) Energy = %.5f +/- %.5f\n", s.mean, en_err); 
    
    double en_exact = SHOExact(path.beta);
    printf("(Exact) Energy = %.5f\n", en_exact);
    
    double err = fabs(s.mean - en_exact) / en_exact;
    printf("Error: %.2f%%\n", err*100.0);
   
     
    if (opt_position_histogram) {
        gsl_histogram *p_histogram;

        size_t nbins = 50;
        p_histogram = gsl_histogram_alloc(nbins);
        gsl_histogram_set_ranges_uniform(p_histogram, -1.0, 1.0);


        for (size_t i = 0; i < trace.positions.count; ++i) {
            gsl_histogram_increment(p_histogram, trace.positions.items[i]);
        }

        double q2_mean = 0.0;
        for (size_t i = 0; i < trace.positions.count; ++i) {
            q2_mean = q2_mean + trace.positions.items[i]*trace.positions.items[i];
        } 
        q2_mean /= trace.positions.count;

        // <q^2> = hbar/(2*m*omega) * coth(lambda/2), where lambda = beta*hbar*omega
        double q2_exact = 0.5/tanh(0.5*path.beta);  

        printf("(PIMC) <q^2> = %.5f\n", q2_mean); 
        printf("(Exact) <q^2> = %.5f\n", q2_exact); 
        err = fabs(q2_mean - q2_exact) / q2_exact;
        printf("Error: %.2f%%\n", err*100.0);

        // draw_histogram("Position histogram", p_histogram);
        // gsl_histogram_free(p_histogram);
    }

    free(trace.energies.items);
    free(trace.positions.items);
    dealloc_beads(&path);    
    close(sockfd);
}

int main(int argc, char *argv[])
{
    char *program_path = shift(&argc, &argv);

    if (argc <= 0) {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: no subcommand is provided\n");
        exit(1);
    }

    const char *subcmd_id = shift(&argc, &argv);
    Subcmd *subcmd = find_subcmd_by_id(subcmds, SUBCMDS_COUNT, subcmd_id);

    if (subcmd != NULL) {
        subcmd->run(program_path, argc, argv);
    } else {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: unknown subcommand  `%s`\n", subcmd_id);
        exit(1);
    }

    arena_free(&arena);

    return 0;
}

int main2()
{
    gsl_histogram *h = gsl_histogram_calloc_uniform(5, 0, 5);
    gsl_histogram_increment(h, 0.5);
    gsl_histogram_increment(h, 1.5);

    gsl_histogram_fprintf(stdout, h, "%f", "%f");

    printf("------------------\n");

    h = gsl_histogram_extend(h);
    gsl_histogram_fprintf(stdout, h, "%f", "%f");

    return 0;

}



// TODO: 
// -- MPI
