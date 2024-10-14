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
} EnergyTrace;

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
    bool verbose = true;
    int sockfd = initServer(verbose);
    if (sockfd < 0) {
        fprintf(stderr, "ERROR: could not start server socket. Exiting...\n");
        exit(1);
    }
    printf("Connection accepted\n");

    bool socket_closed = false; 

    double beta;
    int numTimeSlices;
    int nclients;
    recvDouble(sockfd, &beta);
    recvInt(sockfd, &numTimeSlices);
    recvInt(sockfd, &nclients); 

    double en_exact = SHOExact(beta);
    EnergyTrace tr = {0};

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
    
        Rectangle settingsRect = (Rectangle){ screen_width - screen_width/3, 0, screen_width/3, screen_height };

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
            GuiWindowBox(settingsRect, "Settings");

            int margin = 15;
            Rectangle r = (Rectangle) { settingsRect.x + margin, settingsRect.y + RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT + margin, settingsRect.width, font_size };
            
            GuiLabel(r, TextFormat("clients: %d", nclients));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("beta: %.2f", beta));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("time slices: %d", numTimeSlices));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("Samples: %zu", samples_count));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("Number of bins: %zu", nbins));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("Mean energy: %.5f", mean));
            r.height += font_size + margin;
            
            GuiLabel(r, TextFormat("Exact energy: %.5f", en_exact));
            r.height += font_size + margin;
            
            double exp_error = sigma/sqrt(samples_count);
            GuiLabel(r, TextFormat("Error estimate: %.5f", exp_error));
            r.height += font_size + margin;
            
            double actual_error = mean - en_exact; 
            GuiLabel(r, TextFormat("Actual error: %.5f", actual_error));
            r.height += font_size + margin;
            
            double rel_error = fabs(actual_error) / en_exact;
            GuiLabel(r, TextFormat("Relative error: %.2f%%", rel_error*100.0));
            r.height += font_size + margin;

            /*
            stats_pos.y += font_size; 
            buffer = TextFormat("time slices: %d", numTimeSlices);
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
            */
        }

        EndDrawing();
       
        if (!socket_closed) { 
            // @TODO: recv in non-blocking mode to still have the application responsive
            SocketOpResult r = recvDoubleArray(sockfd, &tr.items, &tr.count);

            if (r == SOCKOP_DISCONNECTED) {
                fprintf(stderr, "Socket closed\n");
                socket_closed = true;
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

