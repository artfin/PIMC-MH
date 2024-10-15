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

typedef struct {
} Props;

#define PROTOCOL_IMPLEMENTATION
#include "protocol.h"

double SHOExact(double beta) {
    //return 0.5/tanh(0.5/T);
    return 0.5/tanh(0.5*beta);
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

int main()
{
    initServer();
    set_verbose(true); 

    int sockfd = 0;
    ConnectionStatus conn = NO_CONNECTION; 
    bool parameters_received = false;

    double beta;
    int numTimeSlices;
    int nclients;
    double en_exact = 0.0;
    
    Trace tr = {0};

    size_t nbins = 20;
    gsl_histogram *h = gsl_histogram_alloc(nbins);

    gsl_histogram_set_ranges_uniform(h, -0.5, 0.5);
    size_t samples_count = 0;

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
    bool editValueBox[2] = { 0 };
    char valTextBox[2][20] = { 0 };

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(COLOR_BACKGROUND);
        
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();
    
        double mean = gsl_histogram_mean(h);
        double sigma = gsl_histogram_sigma(h);

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

        if ((conn == CONNECTION_ESTABLISHED) && !parameters_received) {
            SocketOpResult r;

            r = recvDouble(sockfd, &beta);
            assert(r == SOCKOP_SUCCESS);
            printf("beta received\n");

            r = recvInt(sockfd, &numTimeSlices);
            assert(r == SOCKOP_SUCCESS);
            printf("numTimeSlices received\n");
            
            r = recvInt(sockfd, &nclients);
            assert(r == SOCKOP_SUCCESS);
            printf("nclients received\n");
    
            en_exact = SHOExact(beta);
            parameters_received = true;
            // TODO: what if the parameters won't be recevied in the same tick?
            // just 'assert' that they do.. otherwise we will need to develop some marking scheme 
            // to store the information which parameters have been received and which are still to be awaited 
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
            Vector2 mean_st = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y,
            };
            Vector2 mean_end = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y + world.height,
            };
            DrawLineEx(mean_st, mean_end, 2.0, RED);

            const char *buffer = TextFormat("%.3e", mean);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + (mean - xmin)/(xmax - xmin)*world.width - 0.25*text_len.x,
                world.y + world.height + 0.25 * text_len.y,
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

        if (conn == NO_CONNECTION) {
                
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

        } else if ((conn == CONNECTION_ESTABLISHED) && parameters_received) {
            contentRect.height += font_size + 2*margin;
            GuiLabel(contentRect, TextFormat("Connection established with %s", get_client_ip()));
            contentRect.height += font_size + 2*margin;

            GuiLabel(contentRect, TextFormat("client processes: %d", nclients));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("beta: %.2f", beta));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("time slices: %d", numTimeSlices));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Samples: %zu", samples_count));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Number of bins: %zu", nbins));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Mean energy: %.5f", mean));
            contentRect.height += font_size + margin;
            
            GuiLabel(contentRect, TextFormat("Exact energy: %.5f", en_exact));
            contentRect.height += font_size + margin;
            
            double exp_error = sigma/sqrt(samples_count);
            GuiLabel(contentRect, TextFormat("Error estimate: %.5f", exp_error));
            contentRect.height += font_size + margin;
            
            double actual_error = mean - en_exact; 
            GuiLabel(contentRect, TextFormat("Actual error: %.5f", actual_error));
            contentRect.height += font_size + margin;
            
            double rel_error = fabs(actual_error) / en_exact;
            GuiLabel(contentRect, TextFormat("Relative error: %.2f%%", rel_error*100.0));
            contentRect.height += font_size + margin;
        }

        EndDrawing();
       
        if (conn == CONNECTION_ESTABLISHED) { 
            SocketOpResult r = recvDoubleArray(sockfd, &tr.items, &tr.count);

            if (r == SOCKOP_DISCONNECTED) {
                fprintf(stderr, "Socket closed\n");
                conn = NO_CONNECTION; 
            }
   
            printf("extending histogram with %zu elements\n", tr.count); 
            for (size_t i = 0; i < tr.count; ++i) {
                while ((tr.items[i] < h->range[0] || (tr.items[i] > h->range[h->n]))) {
                    h = gsl_histogram_extend(h);
                    printf("extending histogram: %.5f -- %.5f\n", h->range[0], h->range[h->n]);
                } 
                gsl_histogram_increment(h, tr.items[i]);
            }

            samples_count += tr.count;
        }

        // used as a buffer in recv
        arena_reset(&arena); 
    }

    close(sockfd);
    arena_free(&arena);
    gsl_histogram_free(h);

    return 0;
}

