#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_histogram.h>

#define RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT 32 
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "styles/dark/style_dark.h"

#define COLOR_BACKGROUND GetColor(0x181818FF)
#define FONT_SIZE_LOAD 160 
Font font = {0};

#define COMMON_IMPLEMENTATION
#include "common.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} Trace;

#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 


#define PROTOCOL_IMPLEMENTATION
#include "protocol.h"

// TODO: Send a map (?) with a reference value with respect to which we display the difference

// TODO: Customize the received values from client: we should receive an array of AVERAGES over N values, not raw data 
//       - averageWindow is now read from window
//       - then we have to send it over to client
// TODO: load the default values for histogram ranges from the configuration file

gsl_histogram* gsl_histogram_extend_left(gsl_histogram* h);
gsl_histogram* gsl_histogram_extend_right(gsl_histogram* h);

int main()
{
    initServer();
    set_verbose(true); 

    int result = 0;
    int sockfd = 0;
    ConnectionStatus conn = NO_CONNECTION; 
    bool parameters_exchanged = false;

    // TODO: pack these into a structure?
    double beta, T;
    int numTimeSlices;
    int nclients;
    double refval;
    int blockSize = 1;
    
    Trace tr = {0};

    size_t nbins = 20;
    gsl_histogram *h = gsl_histogram_alloc(nbins);

    gsl_histogram_set_ranges_uniform(h, -0.5, 0.5);
    size_t samples_count = 0;
    size_t packets_count = 0;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "server");
    SetExitKey(KEY_Q);
   
    load_resources();
    int font_size = 24; 

    GuiLoadStyleDark();

    font.baseSize = 100; // @hack: increase the font size for GuiLabel
    GuiSetFont(font);
    font.baseSize = FONT_SIZE_LOAD;  

    bool ylogscale = false;  
    // a pile of values... histogram range & average window size 
    bool editValueBox[3] = { 0 };
    char valTextBox[3][20] = { 0 };

    SetTargetFPS(60);

    double mean = 0.0;
    double var = 0.0;
    double stdev = 0.0;

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(COLOR_BACKGROUND);
        
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();
    
        nbins = h->n;
        double xmin = h->range[0];
        double xmax = h->range[nbins];
        double ymax = gsl_histogram_max_val(h);
        if (ylogscale) ymax = log(ymax);

        int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));
    
        Rectangle settingsRect = (Rectangle){ screen_width - screen_width/3, 0, screen_width/3, screen_height };

        Rectangle world = {
            .x = 0.1*screen_width, 
            .y = 0.1*screen_height, 
            .width = simul_sz, 
            .height = simul_sz
        };
            
        double rect_width = world.width / nbins; 
   
        if (conn == NO_CONNECTION) {
            conn = acceptClientConnection(&sockfd);

            if (conn == CONNECTION_ESTABLISHED) {
                printf("Connection accepted\n");
            }
        } 

        if ((conn == CONNECTION_ESTABLISHED) && !parameters_exchanged) {
            SocketOpResult r;

            r = recvFloat64(sockfd, &beta);
            assert_sockop_result(r); 

            r = recvInt32(sockfd, &numTimeSlices);
            assert_sockop_result(r); 
            
            r = recvInt32(sockfd, &nclients);
            assert_sockop_result(r); 
    
            r = recvFloat64(sockfd, &refval); 
            assert_sockop_result(r); 
                
            sendInt32(sockfd, blockSize);

            parameters_exchanged = true;
  
            // The parameters are exchanged in the BLOCKING mode of the socket
            // and only then we set the socket in the non-blocking mode
                      
            int flags = fcntl(sockfd, F_GETFL, 0);
            flags |= O_NONBLOCK;
            fcntl(sockfd, F_SETFL, flags); 
        } 
        
        for (size_t i = 0; i < nbins; ++i) 
        {
            double factor = 0.9;

            double height;
            if (ylogscale) { 
                height = log(gsl_histogram_get(h, i))/ymax * world.height;
            } else {
                height = gsl_histogram_get(h, i)/ymax * world.height;
            }

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
            DrawLineEx(
                CLITERAL(Vector2){
                    .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                    .y = world.y,
                },
                CLITERAL(Vector2){
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y + world.height,
                }, 2.0, RED);

            const char *buffer = TextFormat("%.3e", mean);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + (mean - xmin)/(xmax - xmin)*world.width - 0.25*text_len.x,
                world.y + world.height + 0.75 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, RED);
        }
        
        {
            const char *buffer = TextFormat("%.3e", xmin);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
        {
            const char *buffer = TextFormat("%.3e", xmax);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + world.width - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
       
        GuiWindowBox(settingsRect, "Settings");
        
        int margin = 15;
        Rectangle contentRect = (Rectangle) { settingsRect.x + margin, settingsRect.y + RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT + margin, settingsRect.width, font_size };
        
        if (conn == NO_CONNECTION) {
            GuiLabel(contentRect, TextFormat("Waiting for client connection..."));
            contentRect.y += font_size;
        }

        GuiCheckBox((Rectangle){ contentRect.x, contentRect.y + contentRect.height, 1.5f*font_size, 1.5f*font_size }, "Y log scale", &ylogscale);
        contentRect.y += font_size + margin;

        if (conn == NO_CONNECTION) 
        {
            GuiLabel((Rectangle){ contentRect.x, contentRect.y + contentRect.height, contentRect.width, font_size }, "Histogram ranges");
            contentRect.y += font_size + margin;
           
            // left range 
            if (GuiTextBox((Rectangle){ contentRect.x, contentRect.y + contentRect.height, 0.4*contentRect.width, 1.5f*font_size }, valTextBox[0], 20, editValueBox[0]))
            {
                editValueBox[0] = !editValueBox[0];

                // Input ended
                if (!editValueBox[0]) {
                    // Try to convert text to float and assign it to the point
                    char *endPtr = NULL;
                    double value = strtod(valTextBox[0], &endPtr);
                    if (endPtr != valTextBox[0]) gsl_histogram_set_ranges_uniform(h, value, xmax);
                }
            }

            // right range 
            if (GuiTextBox((Rectangle){ contentRect.x + contentRect.width/2, contentRect.y + contentRect.height, 0.4*contentRect.width, 1.5f*font_size }, valTextBox[1], 20, editValueBox[1]))
            {
                editValueBox[1] = !editValueBox[1];

                // Input ended
                if (!editValueBox[1])
                {
                    // Try to convert text to float and assign it to the point
                    char *endPtr = NULL;
                    double value = strtod(valTextBox[1], &endPtr);
                    if (endPtr != valTextBox[1]) gsl_histogram_set_ranges_uniform(h, xmin, value);
                }
            }

            // contentRect.y += font_size + margin;
            // GuiLabel((Rectangle){ contentRect.x, contentRect.y + contentRect.height, contentRect.width, font_size }, "Block size (on client)");
            // contentRect.y += font_size + margin;
            // 
            // if (GuiTextBox((Rectangle){ contentRect.x, contentRect.y + contentRect.height, 0.4*contentRect.width, 1.5f*font_size }, valTextBox[2], 20, editValueBox[2]))
            // {
            //     editValueBox[2] = !editValueBox[2];

            //     // Input ended
            //     if (!editValueBox[2]) {
            //         int value = atoi(valTextBox[2]);
            //         if (value != 0) blockSize = value;  
            //     }
            // }

        } else if (parameters_exchanged) {
            contentRect.height += font_size + 2*margin;
            
            if (conn == CONNECTION_ESTABLISHED) {
                GuiLabel(contentRect, TextFormat("Connection established with %s", get_client_ip()));
                contentRect.height += font_size + 2*margin;

                GuiLabel(contentRect, TextFormat("client processes: %d", nclients));
                contentRect.height += font_size + margin;
            }

            T = 1.0 / (beta * Boltzmann_Hartree);
            GuiLabel(contentRect, TextFormat("T: %.2f", T));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("time slices: %d", numTimeSlices));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Samples: %zu", samples_count));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Packets: %zu", packets_count));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Number of bins: %zu", nbins));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Mean: %.5e", mean));
            contentRect.height += font_size + margin;
            
            double exp_error = stdev/sqrt(samples_count);
            GuiLabel(contentRect, TextFormat("Error estimate: %.5e", exp_error));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Reference: %.5e", refval));
            contentRect.height += font_size + margin;
            
            double actual_error = mean - refval; 
            GuiLabel(contentRect, TextFormat("Actual error: %.5e", actual_error));
            contentRect.height += font_size + margin;
            
            double rel_error = fabs(actual_error) / refval;
            GuiLabel(contentRect, TextFormat("Relative error: %.3f%%", rel_error*100.0));
            contentRect.height += font_size + margin;
        }

        EndDrawing();
       
        if (conn == CONNECTION_ESTABLISHED) { 
            SocketOpResult r = recvFloat64Array(sockfd, &tr.items, &tr.count);

            if (r == SOCKOP_DISCONNECTED) {
                fprintf(stderr, "Socket closed\n");
                conn = DISCONNECTED; 
            }

            if (r == SOCKOP_SUCCESS) { // otherwise we are waiting for the data packet to arrive
                printf("extending histogram with %zu elements\n", tr.count); 

                double packet_mean = 0.0;
                for (size_t i = 0; i < tr.count; ++i) {
                    while (tr.items[i] < h->range[0]) {
                        h = gsl_histogram_extend_left(h);
                        printf("extending histogram: %.5f -- %.5f\n", h->range[0], h->range[h->n]);
                    }

                    while (tr.items[i] > h->range[h->n]) {
                        h = gsl_histogram_extend_right(h);
                        printf("extending histogram: %.5f -- %.5f\n", h->range[0], h->range[h->n]);
                    }

                    gsl_histogram_increment(h, tr.items[i]);
                    packet_mean += tr.items[i];
                }
                packet_mean /= tr.count;

                double packet_var  = 0.0;
                for (size_t i = 0; i < tr.count; ++i) {
                    packet_var += (tr.items[i] - packet_mean)*(tr.items[i] - packet_mean);
                }
                packet_var /= (tr.count - 1);

                assert(blockSize == 1); 
                var = ((samples_count-1)*var + (tr.count-1)*packet_var) / (samples_count+tr.count-1) + samples_count*tr.count*(mean - packet_mean)*(mean - packet_mean)/(samples_count+tr.count)/(samples_count+tr.count - 1);
                mean = (mean*samples_count + tr.count*packet_mean) / (samples_count + tr.count);
                stdev = sqrt(var);

                samples_count += tr.count;
                packets_count++;
            }
        }

        // used as a buffer in recv
        arena_reset(&arena); 
    }

defer:
    if (sockfd) {
        printf("Closing socket...\n");
        close(sockfd);
    }

    arena_free(&arena);
    gsl_histogram_free(h);

    return result;
}

gsl_histogram* gsl_histogram_extend_left(gsl_histogram* h)
{
    size_t nbins = h->n;
    double add_bins = 1;
    
    double xmin = h->range[0];
    double xmax = h->range[nbins];
    double dx = h->range[1] - h->range[0];
    
    double new_xmin = xmin - add_bins*dx; 
    
    gsl_histogram *new_h = gsl_histogram_alloc(nbins + add_bins);
    gsl_histogram_set_ranges_uniform(new_h, new_xmin, xmax);
        
    size_t nc = add_bins; // cursor over the new histogram
    for (size_t i = 0; i < nbins; ++i) {
        new_h->bin[nc++] = gsl_histogram_get(h, i); 
    }

    gsl_histogram_free(h);

    return new_h;
}

gsl_histogram* gsl_histogram_extend_right(gsl_histogram* h)
{
    size_t nbins = h->n;
    double add_bins = 1;
    
    double xmin = h->range[0];
    double xmax = h->range[nbins];
    double dx = h->range[1] - h->range[0];
    
    double new_xmax = xmax + add_bins*dx; 
    
    gsl_histogram *new_h = gsl_histogram_alloc(nbins + add_bins);
    gsl_histogram_set_ranges_uniform(new_h, xmin, new_xmax);
        
    size_t nc = 0; // cursor over the new histogram
    for (size_t i = 0; i < nbins; ++i) {
        new_h->bin[nc++] = gsl_histogram_get(h, i); 
    }

    gsl_histogram_free(h);

    return new_h;
}

