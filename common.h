#ifndef COMMON_H_
#define COMMON_H_

#define return_defer(value) do { result = (value); goto defer; } while (0)

#define DA_INIT_CAP 256

#define da_append(da, item)                                                          \
    do {                                                                             \
        if ((da)->count >= (da)->capacity) {                                         \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2;   \
            (da)->items = realloc((da)->items, (da)->capacity*sizeof(*(da)->items)); \
            assert((da)->items != NULL && "Buy more RAM lol");                       \
        }                                                                            \
                                                                                     \
        (da)->items[(da)->count++] = (item);                                         \
    } while (0)

#define DEFINE_SUBCMD(name, desc) \
    {                             \
        .run = subcmd_##name,     \
        .id  = #name,             \
        .description = desc       \
    } 

typedef struct {
    void (*run)(const char *program_path, int argc, char **argv);
    const char *id;
    const char *description;
} Subcmd;

typedef struct {
    int rank;
    int size;
} MPI_Context;

typedef struct {
    void (*run)(MPI_Context ctx, int argc, char **argv);
    const char *id;
    const char *description;
} MPI_Subcmd;


Subcmd *find_subcmd_by_id(Subcmd *subcmds, size_t subcmds_count, const char *id); 
void usage(const char *program_path, Subcmd *subcmds, size_t subcmds_count);
char* shift(int *argc, char ***argv);

typedef struct {
    size_t CenterOfMass;
    size_t Staging;
} AcceptanceRate;

typedef struct {
    double **beads;
    size_t numParticles;
    size_t numTimeSlices;
    double tau;
    double beta;
} Path;

#ifdef COMMON_IMPLEMENTATION

Subcmd *find_subcmd_by_id(Subcmd *subcmds, size_t subcmds_count, const char *id) 
{
    for (size_t k = 0; k < subcmds_count; ++k) {
        if (strcmp(subcmds[k].id, id) == 0) {
            return &subcmds[k];
        }
    }

    return NULL; 
}

void usage(const char *program_path, Subcmd *subcmds, size_t subcmds_count)
{
    fprintf(stderr, "Usage: %s [subcommand]\n", program_path);
    fprintf(stderr, "Subcommands:\n");

    int width = 0;
    for (size_t k = 0; k < subcmds_count; ++k) {
        int len = strlen(subcmds[k].id);
        if (width < len) width = len;
    }

    for (size_t k = 0; k < subcmds_count; ++k) {
        fprintf(stderr, "    %-*s - %s\n", width, subcmds[k].id, subcmds[k].description);
    }
}

char* shift(int *argc, char ***argv)
{
    assert(*argc > 0);
    char *result = *argv[0];

    *argc -= 1;
    *argv += 1;

    return result; 
}


#endif // COMMON_IMPLEMENTATION

#endif // COMMON_H_
