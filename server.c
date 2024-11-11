#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <float.h>

#include <gsl/gsl_histogram.h>

#define RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT 32 
#define RAYGUI_IMPLEMENTATION
// #define RAYGUI_NO_ICONS
#include "raygui.h"
#include "styles/dark/style_dark.h"

#define COLOR_BACKGROUND GetColor(0x181818FF)
#define FONT_SIZE_LOAD 160 
Font font = {0};
int FONT_SIZE = 24; 

#define COMMON_IMPLEMENTATION
#include "common.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0}; 
static Arena arena_str = {0};

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} Trace;

typedef struct {
    double min;
    double max;
} Range; 

#define PROTOCOL_IMPLEMENTATION
#include "protocol.h"

typedef struct {
    char name[MAX_NAME_SIZE];
    
    gsl_histogram *h;
    size_t samples_count; 
    
    // NOTE: 
    // We have three sets of minimum/maximum values:
    // 1) displayed region (display_xmin, display_xmax)
    // 2) stored in histogram (hist_xmin, hist_xmax)
    // 3) actual data (data_xmin, data_xmax)
    Range data_range;
    Range display_range;
    Range hist_range;
  
    double refval;  
    double mean;
    double var;
    double stdev;

    size_t packets_count;
} Tab;

#define MAX_TABS 3
static char tab_names[MAX_TABS][MAX_NAME_SIZE] = {0};
static char **ptab_names = NULL; // different representation of tab_names for GuiTabBar
static Tab TABS[MAX_TABS] = {0};
static size_t TAB_COUNT = 0;

#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 


#define return_defer(value) do { result = (value); goto defer; } while (0)

// TODO: Send a map with a reference value with respect to which we display the difference

// TODO: Customize the received values from client: we should receive an array of AVERAGES over N values, not raw data 
//       - averageWindow is now read from window
//       - then we have to send it over to client
// TODO: load the default values for histogram ranges from the configuration file

gsl_histogram* gsl_histogram_extend_left(gsl_histogram* h);
gsl_histogram* gsl_histogram_extend_right(gsl_histogram* h);
void set_histogram_ranges_from_packet(gsl_histogram *h, Trace packet);

bool get_histogram_ranges(Rectangle contentRect, const char *description, double *lhs, double *rhs, bool editMode)
{
    bool result = false;
    
    static bool editValueBox[2] = {0};
    static char valTextBox[2][20] = {0};
    
    int margin = 15;
    
    if (!editMode) {
        guiState = STATE_DISABLED;
    }

    GuiLabel((Rectangle){ contentRect.x, contentRect.y + contentRect.height, contentRect.width, FONT_SIZE }, description);
    contentRect.y += FONT_SIZE + margin;

    // left range 
    if (GuiTextBox((Rectangle){ contentRect.x, contentRect.y + contentRect.height, 0.4*contentRect.width, 1.5f*FONT_SIZE }, valTextBox[0], 20, editValueBox[0]))
    {
        editValueBox[0] = !editValueBox[0];

        if (!editValueBox[0]) {
            // Try to convert text to float and assign it to the point
            char *endPtr = NULL;
            double value = strtod(valTextBox[0], &endPtr);
            if (endPtr != valTextBox[0]) *lhs = value;
            return_defer(true);
        }
    }

    // right range 
    if (GuiTextBox((Rectangle){ contentRect.x + contentRect.width/2, contentRect.y + contentRect.height, 0.4*contentRect.width, 1.5f*FONT_SIZE }, valTextBox[1], 20, editValueBox[1]))
    {
        editValueBox[1] = !editValueBox[1];

        if (!editValueBox[1]) {
            // Try to convert text to float and assign it to the point
            char *endPtr = NULL;
            double value = strtod(valTextBox[1], &endPtr);
            if (endPtr != valTextBox[1]) *rhs = value;
            return_defer(true);
        }
    }

defer:
    guiState = STATE_NORMAL;
    return result; 
}

void display_histogram(Rectangle r, gsl_histogram *h, Range display_range, bool ylogscale)
{
    size_t nbins = h->n;
    assert(nbins > 1);
    double xmin = h->range[0];
    double dx = h->range[1] - h->range[0];
    
    double ymax = FLT_MIN; // this the maximum "y" in the displayed region
    size_t displayed_nbins = 0;

    for (size_t i = 0; i < nbins; ++i) {
        double x = xmin + i*dx;
        if (x < display_range.min) continue;
        if (x + dx > display_range.max) break;

        double y = gsl_histogram_get(h, i);
        if (y > ymax) ymax = y;

        displayed_nbins++;
    }

    if (ylogscale) ymax = log(ymax);

    double col_width = (displayed_nbins > 0) ? r.width / displayed_nbins : 0.0; 

    size_t display_index = 0;
    for (size_t i = 0; i < nbins; ++i) 
    {
        double x = xmin + i*dx;
        if (x < display_range.min) continue;
        if (x + dx > display_range.max) break;

        double scale_height = 0.9;

        double height;
        if (ylogscale) { 
            height = log(1.0f + gsl_histogram_get(h, i))/ymax * r.height;
        } else {
            height = gsl_histogram_get(h, i)/ymax * r.height;
        }

        Rectangle col = {
            .x = r.x + display_index*col_width,
            .y = r.y + r.height - scale_height*height, 
            .width = col_width,
            .height = scale_height*height, 
        };

        DrawRectangleLinesEx(col, 2.0, YELLOW);
        display_index++;
    }
}

void display_mean(Rectangle world, double mean, Range display_range)
{
    DrawLineEx(
        CLITERAL(Vector2){
        .x = world.x + (mean - display_range.min)/(display_range.max - display_range.min)*world.width, 
        .y = world.y,
        },
        CLITERAL(Vector2){
        .x = world.x + (mean - display_range.min)/(display_range.max - display_range.min)*world.width, 
        .y = world.y + world.height,
        }, 2.0, RED);

    const char *buffer = TextFormat("%.3e", mean);
    Vector2 text_len = MeasureTextEx(font, buffer, FONT_SIZE, 0);
    Vector2 text_pos = {
        world.x + (mean - display_range.min)/(display_range.max - display_range.min)*world.width - 0.25*text_len.x,
        world.y + world.height + 0.75 * text_len.y,
    };
    DrawTextEx(font, buffer, text_pos, FONT_SIZE, 0, RED);
} 

void display_xlabels(Rectangle world, Range display_range)
{
    {
        const char *buffer = TextFormat("%.3e", display_range.min);
        Vector2 text_len = MeasureTextEx(font, buffer, FONT_SIZE, 0);
        Vector2 text_pos = {
            world.x - 0.25 * text_len.x,
            world.y + world.height + 0.25 * text_len.y,
        };
        DrawTextEx(font, buffer, text_pos, FONT_SIZE, 0, WHITE);
    }

    {
        const char *buffer = TextFormat("%.3e", display_range.max);
        Vector2 text_len = MeasureTextEx(font, buffer, FONT_SIZE, 0);
        Vector2 text_pos = {
            world.x + world.width - 0.25 * text_len.x,
            world.y + world.height + 0.25 * text_len.y,
        };
        DrawTextEx(font, buffer, text_pos, FONT_SIZE, 0, WHITE);
    }
}

void show_grid(Rectangle world)
{
    size_t nlines = 10;
    double DEFAULT_LW = 0.5;

    for (size_t i = 0; i < nlines; ++i) 
    {
        double lw = (i == (nlines + 1)/2) ? 3*DEFAULT_LW : DEFAULT_LW;
        Color color = (i == (nlines + 1)/2) ? GetColor(0xE1BFF264) : GetColor(0xC8C8C864); 

        double dx = (float)i/nlines*world.width; 
        DrawLineEx(
            CLITERAL(Vector2){
            .x = world.x + dx, 
            .y = world.y,
            },
            CLITERAL(Vector2){
            .x = world.x + dx, 
            .y = world.y + world.height,
            }, lw, color);
        
        double dy = (float)i/nlines*world.height; 
        DrawLineEx(
            CLITERAL(Vector2){
            .x = world.x, 
            .y = world.y + dy,
            },
            CLITERAL(Vector2){
            .x = world.x + world.width, 
            .y = world.y + dy,
            }, lw, color);
    }
}

char **convert_char2d_to_ptr(char (*names)[MAX_NAME_SIZE], size_t size);

Tab* tab_alloc(const char *name) {
    assert(TAB_COUNT < MAX_TABS);
   
    printf("INFO: Creating tab (%zu) with name = %s\n", TAB_COUNT, name);
    Tab *tab = &TABS[TAB_COUNT];

    assert(strlen(name) < MAX_NAME_SIZE); 
    strcpy(tab->name, name);
    strcpy(tab_names[TAB_COUNT], name);

    size_t nbins = 20;
    tab->h = gsl_histogram_alloc(nbins);
    gsl_histogram_set_ranges_uniform(tab->h, -0.5, 0.5);
    
    tab->data_range = (Range) { 
        .min = FLT_MAX, 
        .max = FLT_MIN 
    };
    
    tab->display_range = (Range) { 
        .min = tab->h->range[0],
        .max = tab->h->range[nbins]
    };

    tab->hist_range = (Range) {
        .min = tab->h->range[0],
        .max = tab->h->range[nbins],
    };

    TAB_COUNT++;

    ptab_names = convert_char2d_to_ptr(tab_names, TAB_COUNT);

    return tab;
}


void add_packet_to_tab(Tab *tab, Trace tr)
{
    printf("INFO: Extending histogram %s with %zu elements\n", tab->name, tr.count); 

    if (tab->packets_count == 0) {
        set_histogram_ranges_from_packet(tab->h, tr);
        printf("INFO: Setting initial range for histogram %s to [%.5lf...%.5lf]\n", 
                tab->name, tab->h->range[0], tab->h->range[tab->h->n]); 
    }

    double packet_mean = 0.0;
    for (size_t i = 0; i < tr.count; ++i) {
        while (tr.items[i] < tab->h->range[0]) {
            tab->h = gsl_histogram_extend_left(tab->h);
            printf("extending histogram: %.5f -- %.5f\n", tab->h->range[0], tab->h->range[tab->h->n]);
        }

        while (tr.items[i] > tab->h->range[tab->h->n]) {
            tab->h = gsl_histogram_extend_right(tab->h);
            printf("extending histogram: %.5f -- %.5f\n", tab->h->range[0], tab->h->range[tab->h->n]);
        }

        gsl_histogram_increment(tab->h, tr.items[i]);
        packet_mean += tr.items[i];

        if (tr.items[i] > tab->data_range.max) tab->data_range.max = tr.items[i];
        if (tr.items[i] < tab->data_range.min) tab->data_range.min = tr.items[i];
    }
    packet_mean /= tr.count;

    double packet_var  = 0.0;
    for (size_t i = 0; i < tr.count; ++i) {
        packet_var += (tr.items[i] - packet_mean)*(tr.items[i] - packet_mean);
    }
    packet_var /= (tr.count - 1);

    tab->var = ((tab->samples_count-1)*tab->var + (tr.count-1)*packet_var) / (tab->samples_count + tr.count-1) + 
                 tab->samples_count*tr.count*(tab->mean - packet_mean)*(tab->mean - packet_mean)/(tab->samples_count + tr.count)/(tab->samples_count + tr.count - 1);
    tab->mean = (tab->mean*tab->samples_count + tr.count*packet_mean) / (tab->samples_count + tr.count);
    tab->stdev = sqrt(tab->var);

    tab->samples_count += tr.count;
    tab->packets_count++;
}


void test(const char** text, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        printf("text[%zu] = %s\n", i, text[i]);
    }
}

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
    
    Trace tr = {0};
    // Trace necklace_sizes = {0};

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "server");
    SetExitKey(KEY_Q);
   
    load_resources();

    GuiLoadStyleDark();
    GuiSetStyle(CHECKBOX, TEXT_PADDING, 12);

    font.baseSize = 100; // @hack: increase the font size for GuiLabel
    GuiSetFont(font);
    font.baseSize = FONT_SIZE_LOAD;  

    bool ylogscale = false;  
    bool xadaptive = true;

    SetTargetFPS(60);

    strcpy(tab_names[0], "Necklace size");
    strcpy(tab_names[1], "Energy");
    for (size_t i = 0; i < MAX_TABS; ++i) {
        char *name = tab_names[i];
        if (strlen(name) > 0) {
            tab_alloc(name);
        }
    }
     
    Tab *active_tab = &TABS[0]; 

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(COLOR_BACKGROUND);
        
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();

        active_tab->hist_range.min = active_tab->h->range[0];
        active_tab->hist_range.max = active_tab->h->range[active_tab->h->n]; 
        
        if (xadaptive) {
            active_tab->display_range.min = active_tab->hist_range.min;
            active_tab->display_range.max = active_tab->hist_range.max;
        }

        int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));
    
        Rectangle settingsRect = (Rectangle){ screen_width - screen_width/3, 0, screen_width/3, screen_height };

        Rectangle world = {
            .x = 0.1*screen_width, 
            .y = 0.1*screen_height, 
            .width = simul_sz, 
            .height = simul_sz
        };
        DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);

        Rectangle bar = CLITERAL(Rectangle) {
            .x = 0.02*screen_width, 
            .y = 0.02*screen_height, 
            .width = simul_sz, 
            .height = 50 
        };
   
        // GuiSetStyle(TOGGLE, TEXT_COLOR_NORMAL, 0xff0000ff);
        int cursor;
        GuiTabBar(bar, (const char**) ptab_names, TAB_COUNT, &cursor);
        active_tab = &TABS[cursor]; 
   
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
  
            // TODO: can we make this more flexible? 
            // can we communicate the reference value with the data packet?
            for (size_t i = 0; i < TAB_COUNT; ++i) {
                Tab *tab = &TABS[i];

                char *name = NULL;
                double refval;
                r = recvNamedFloat64(sockfd, &name, &refval);
                assert_sockop_result(r); 
                
                if (strcmp(tab->name, name) != 0) {
                    printf("ERROR: Expected '%s' but got '%s'\n", tab->name, name);
                    return_defer(1); 
                }
                tab->refval = refval;
            }
                
            parameters_exchanged = true;
  
            // The parameters are exchanged in the BLOCKING mode of the socket
            // and only then we set the socket in the non-blocking mode
            int flags = fcntl(sockfd, F_GETFL, 0);
            flags |= O_NONBLOCK;
            fcntl(sockfd, F_SETFL, flags); 
        } 
        
        show_grid(world);
        display_histogram(world, active_tab->h, active_tab->display_range, ylogscale);

        if (active_tab->samples_count > 0) { 
            display_mean(world, active_tab->mean, active_tab->display_range);      
        }

        display_xlabels(world, active_tab->display_range);

       
        GuiWindowBox(settingsRect, "Settings");
        
        int margin = 15;
        Rectangle contentRect = (Rectangle) { settingsRect.x + margin, settingsRect.y + RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT + margin, settingsRect.width, FONT_SIZE };
        
        if (conn == NO_CONNECTION) {
            GuiLabel(contentRect, TextFormat("Waiting for client connection..."));
            contentRect.y += FONT_SIZE;
        } else if (conn == CONNECTION_ESTABLISHED) {
            GuiLabel(contentRect, TextFormat("Connection established with %s", get_client_ip()));
            contentRect.y += FONT_SIZE;
        }

        {
            GuiCheckBox((Rectangle){ contentRect.x, contentRect.y + contentRect.height, 1.5f*FONT_SIZE, 1.5f*FONT_SIZE }, "Y logscale", &ylogscale);
            GuiCheckBox((Rectangle){ contentRect.x + 0.5*contentRect.width, contentRect.y + contentRect.height, 1.5f*FONT_SIZE, 1.5f*FONT_SIZE }, "X adaptive", &xadaptive);
            contentRect.y += FONT_SIZE + margin;
        }


        // TODO: display boxes for ranges of histogram all the time
        // we could set INITIAL ranges and then modify the displayed range of the histogram
        if (conn == NO_CONNECTION) {
            double lhs = active_tab->hist_range.min;
            double rhs = active_tab->hist_range.max;

            if (get_histogram_ranges(contentRect, "Initial range:", &lhs, &rhs, true)) {
                if (lhs < rhs) {
                    gsl_histogram_set_ranges_uniform(active_tab->h, lhs, rhs);
                    active_tab->display_range.min = lhs;
                    active_tab->display_range.max = rhs; 

                    GuiSetStyle(TEXTBOX, BORDER_COLOR_PRESSED, 0x000000ff);
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_FOCUSED, 0xe1e1e1ff);  
                } else {
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_PRESSED, 0xff0000ff);
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_FOCUSED, 0xff0000ff);  
                }
            }
        } else if (conn == CONNECTION_ESTABLISHED) {
            double lhs = active_tab->display_range.min;
            double rhs = active_tab->display_range.max;

            if (get_histogram_ranges(contentRect, "Display range:", &lhs, &rhs, !xadaptive)) {
                if (lhs < rhs) {
                    active_tab->display_range.min = lhs;
                    active_tab->display_range.max = rhs; 
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_PRESSED, 0x000000ff);
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_FOCUSED, 0xe1e1e1ff); 
                } else {
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_PRESSED, 0xff0000ff);
                    GuiSetStyle(TEXTBOX, BORDER_COLOR_FOCUSED, 0xff0000ff);  
                }
            }
        }

        contentRect.y += contentRect.height; 

        if (parameters_exchanged) {
            contentRect.y += FONT_SIZE + 4*margin;
            
            GuiLabel(contentRect, TextFormat("client processes: %d", nclients));
            contentRect.y += 0.9*FONT_SIZE;

            T = 1.0 / (beta * Boltzmann_Hartree);
            GuiLabel(contentRect, TextFormat("T: %.2f", T));
            contentRect.y += 0.9*FONT_SIZE;
            
            GuiLabel(contentRect, TextFormat("time slices: %d", numTimeSlices));
            contentRect.y += 0.9*FONT_SIZE;
            
            GuiLabel(contentRect, TextFormat("Samples: %zu", active_tab->samples_count));
            contentRect.y += 0.9*FONT_SIZE;
            
            GuiLabel(contentRect, TextFormat("Packets: %zu", active_tab->packets_count));
            contentRect.y += 0.9*FONT_SIZE;
           
            GuiLabel(contentRect, TextFormat("Data min: %.3e", active_tab->data_range.min));
            contentRect.y += 0.9*FONT_SIZE;

            GuiLabel(contentRect, TextFormat("Data max: %.3e", active_tab->data_range.max));
            contentRect.y += 0.9*FONT_SIZE;

            GuiLabel(contentRect, TextFormat("Number of bins: %zu", active_tab->h->n));
            contentRect.y += 0.9*FONT_SIZE;
            
            GuiLabel(contentRect, TextFormat("Mean: %.5e", active_tab->mean));
            contentRect.y += 0.9*FONT_SIZE;
            
            double exp_error = active_tab->stdev/sqrt(active_tab->samples_count);
            GuiLabel(contentRect, TextFormat("Error estimate: %.5e", exp_error));
            contentRect.y += 0.9*FONT_SIZE;
            
            GuiLabel(contentRect, TextFormat("Reference: %.5e", active_tab->refval));
            contentRect.y += 0.9*FONT_SIZE;
            
            double actual_error = active_tab->mean - active_tab->refval; 
            GuiLabel(contentRect, TextFormat("Actual error: %.5e", actual_error));
            contentRect.y += 0.9*FONT_SIZE;
            
            double rel_error = fabs(actual_error) / active_tab->refval;
            GuiLabel(contentRect, TextFormat("Relative error: %.3f%%", rel_error*100.0));
            contentRect.y += 0.9*FONT_SIZE;
        }

        EndDrawing();
      
        if (conn == CONNECTION_ESTABLISHED) {
            SocketOpResult r;

            char *packet_name = NULL;
            r = recvNamedFloat64Array(sockfd, &packet_name, &tr.items, &tr.count);

            if (r == SOCKOP_DISCONNECTED) {
                fprintf(stderr, "Socket closed\n");
                conn = DISCONNECTED;
                continue; 
            } else if (r == SOCKOP_SUCCESS) { // otherwise we are waiting for the data packet to arrive
                // NOTE: "-1" is added so that strncpy could add a null-terminator..   
                // strncpy((char*) tab_names[0], name, MAX_NAME_SIZE - 1); 
                // printf("name: %s\n", tab_names[0]);

                Tab *tab = NULL;
               
                if ((TAB_COUNT == 1) && (strcmp(active_tab->name, "default") == 0)) {
                    tab = active_tab;
                    strcpy(tab->name, packet_name);
                    strcpy(tab_names[0], packet_name);
                } else {
                    for (size_t i = 0; i < MAX_TABS; ++i) {
                        if (strcmp(TABS[i].name, packet_name) == 0) {
                            tab = &TABS[i];
                        }
                    }

                    if (tab == NULL) {
                        tab = tab_alloc(packet_name);
                    }
                }

                add_packet_to_tab(tab, tr);
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
    arena_free(&arena_str);

    // TODO: free tabs
    //gsl_histogram_free(h);

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

void set_histogram_ranges_from_packet(gsl_histogram *h, Trace packet)
{
    double packet_min = FLT_MAX;
    double packet_max = FLT_MIN;

    for (size_t i = 0; i < packet.count; ++i) {
        if (packet.items[i] < packet_min) packet_min = packet.items[i];
        if (packet.items[i] > packet_max) packet_max = packet.items[i];
    } 

    gsl_histogram_set_ranges_uniform(h, packet_min, packet_max);
}

char **convert_char2d_to_ptr(char (*names)[MAX_NAME_SIZE], size_t size) {
    char **p = (char **) arena_alloc(&arena_str, size * sizeof(char *));
    for (size_t i = 0; i < size; i++) {
        p[i] = names[i];
    }

    return p;
}

