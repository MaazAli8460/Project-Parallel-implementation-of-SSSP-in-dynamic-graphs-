#define _POSIX_C_SOURCE 200112L
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <mpi.h>
#include <metis.h>
#include <omp.h>

#define INF INT_MAX
#define BATCH_SIZE 30 

typedef struct { int to, w; } Edge;
typedef struct {
    int V;
    Edge** adj;
    int* adj_sizes;
    int* adj_caps;
    char* eSet;
    int eSet_size;
} Graph;


typedef struct {
    int* dist;
    int* parent;
    char* aff;
    char* affDel;
    int size;
} SSSP;


typedef struct {
    double comm_time;   
    double cpu_time;     
    double total_cpu_time; 
    double wall_time;   
    double del_time;     
    double ins_time;     
    double prop_time;    
    double relax_time;   
} Timing;

/*---------------- Graph routines ----------------*/
Graph* initGraph(int n) {
    Graph* G = malloc(sizeof(Graph));
    if (!G) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    G->V = n;
    G->adj = calloc(n, sizeof(Edge*));
    G->adj_sizes = calloc(n, sizeof(int));
    G->adj_caps = calloc(n, sizeof(int));
    G->eSet_size = 1000003;
    G->eSet = calloc(G->eSet_size, 1);
    if (!G->adj || !G->adj_sizes || !G->adj_caps || !G->eSet) {
        fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    return G;
}

void freeGraph(Graph* G) {
    for (int i = 0; i < G->V; i++) free(G->adj[i]);
    free(G->adj); free(G->adj_sizes); free(G->adj_caps); free(G->eSet); free(G);
}

void resizeGraph(Graph* G, int n) {
    if (n <= G->V) return;
    G->adj = realloc(G->adj, n * sizeof(Edge*));
    G->adj_sizes = realloc(G->adj_sizes, n * sizeof(int));
    G->adj_caps = realloc(G->adj_caps, n * sizeof(int));
    if (!G->adj || !G->adj_sizes || !G->adj_caps) {
        fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    for (int i = G->V; i < n; i++) {
        G->adj[i] = NULL;
        G->adj_sizes[i] = 0;
        G->adj_caps[i] = 0;
    }
    G->V = n;
}

unsigned long long edgeKey(int u, int v) {
    int mn = u < v ? u : v;
    int mx = u < v ? v : u;
    return ((unsigned long long)mn << 32) | (unsigned int)mx;
}

void addEdge(Graph* G, int u, int v, int w) {
    if (u >= G->V || v >= G->V) resizeGraph(G, (u > v ? u : v) + 1);
    if (G->adj_sizes[u] >= G->adj_caps[u]) {
        G->adj_caps[u] = G->adj_caps[u] ? G->adj_caps[u] * 2 : 4;
        G->adj[u] = realloc(G->adj[u], G->adj_caps[u] * sizeof(Edge));
        if (!G->adj[u]) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    }
    G->adj[u][G->adj_sizes[u]++] = (Edge){v, w};
    if (G->adj_sizes[v] >= G->adj_caps[v]) {
        G->adj_caps[v] = G->adj_caps[v] ? G->adj_caps[v] * 2 : 4;
        G->adj[v] = realloc(G->adj[v], G->adj_caps[v] * sizeof(Edge));
        if (!G->adj[v]) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    }
    G->adj[v][G->adj_sizes[v]++] = (Edge){u, w};
    G->eSet[edgeKey(u, v) % G->eSet_size] = 1;
}

int edgeExists(Graph* G, int u, int v) {
    return (u < G->V && v < G->V) ? G->eSet[edgeKey(u, v) % G->eSet_size] : 0;
}

void removeEdge(Graph* G, int u, int v) {
    if (u >= G->V || v >= G->V) return;
    for (int i = 0; i < G->adj_sizes[u]; i++) {
        if (G->adj[u][i].to == v) {
            G->adj[u][i] = G->adj[u][--G->adj_sizes[u]];
            break;
        }
    }
    for (int i = 0; i < G->adj_sizes[v]; i++) {
        if (G->adj[v][i].to == u) {
            G->adj[v][i] = G->adj[v][--G->adj_sizes[v]];
            break;
        }
    }
    G->eSet[edgeKey(u, v) % G->eSet_size] = 0;
}

SSSP* initSSSP(int n) {
    SSSP* T = malloc(sizeof(SSSP));
    if (!T) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    T->size = n;
    T->dist = malloc(n * sizeof(int));
    T->parent = malloc(n * sizeof(int));
    T->aff = calloc(n, 1);
    T->affDel = calloc(n, 1);
    if (!T->dist || !T->parent || !T->aff || !T->affDel) {
        fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    for (int i = 0; i < n; i++) {
        T->dist[i] = INF;
        T->parent[i] = -1;
    }
    return T;
}

void freeSSSP(SSSP* T) {
    free(T->dist); free(T->parent); free(T->aff); free(T->affDel); free(T);
}

void resizeSSSP(SSSP* T, int n) {
    if (n <= T->size) return;
    T->dist = realloc(T->dist, n * sizeof(int));
    T->parent = realloc(T->parent, n * sizeof(int));
    T->aff = realloc(T->aff, n);
    T->affDel = realloc(T->affDel, n);
    if (!T->dist || !T->parent || !T->aff || !T->affDel) {
        fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    for (int i = T->size; i < n; i++) {
        T->dist[i] = INF;
        T->parent[i] = -1;
        T->aff[i] = 0;
        T->affDel[i] = 0;
    }
    T->size = n;
}

int bfs_component(const Graph* G, int src, char* reachable) {
    int* queue = malloc(G->V * sizeof(int));
    if (!queue) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    int front = 0, rear = 0, size_reachable = 0;
    
    memset(reachable, 0, G->V);
    reachable[src] = 1;
    queue[rear++] = src;
    size_reachable++;
    
    while (front < rear) {
        int u = queue[front++];
        for (int i = 0; i < G->adj_sizes[u]; i++) {
            int v = G->adj[u][i].to;
            if (!reachable[v]) {
                reachable[v] = 1;
                queue[rear++] = v;
                size_reachable++;
            }
        }
    }
    free(queue);
    return size_reachable;
}


typedef struct { int d, v; } PQE;
static int cmpPq(const void* a, const void* b) {
    return ((PQE*)a)->d - ((PQE*)b)->d;
}

void dijkstra(const Graph* G, SSSP* T, int src) {
    resizeSSSP(T, G->V);
    for (int i = 0; i < G->V; i++) {
        T->dist[i] = INF;
        T->parent[i] = -1;
    }
    T->dist[src] = 0;
    PQE* pq = malloc(G->V * sizeof(PQE));
    if (!pq) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    int psz = 0;
    pq[psz++] = (PQE){0, src};
    while (psz) {
        qsort(pq, psz, sizeof(PQE), cmpPq);
        PQE e = pq[0];
        memmove(pq, pq + 1, (--psz) * sizeof(PQE));
        if (e.d != T->dist[e.v]) continue;
        for (int i = 0; i < G->adj_sizes[e.v]; i++) {
            Edge ed = G->adj[e.v][i];
            int nd = (e.d == INF || ed.w == INF) ? INF : e.d + ed.w;
            if (nd < T->dist[ed.to]) {
                T->dist[ed.to] = nd;
                T->parent[ed.to] = e.v;
                pq[psz++] = (PQE){nd, ed.to};
            }
        }
    }
    free(pq);
}

int processDeletions(Graph* G, SSSP* T, int del[][2], int delSz, int* affCount) {
    int localAff = 0;
    #pragma omp parallel for schedule(static) reduction(+:localAff)
    for (int i = 0; i < delSz; i++) {
        int u = del[i][0], v = del[i][1];
        if (u >= G->V || v >= G->V) continue;
        int inTree = (T->parent[u] == v) || (T->parent[v] == u);
        if (inTree) {
            int c = (T->parent[v] == u ? v : u);
            T->dist[c] = INF;
            T->parent[c] = -1;
            T->aff[c] = 1;
            T->affDel[c] = 1;
            localAff++;
        }
        if (edgeExists(G, u, v)) removeEdge(G, u, v);
    }
    *affCount += localAff;
    return localAff > 0;
}

int processInsertions(Graph* G, SSSP* T, int ins[][2], int insSz, int* wgt, int* affCount) {
    int mx = G->V;
    for (int i = 0; i < insSz; i++) {
        mx = (ins[i][0] > mx ? ins[i][0] : mx);
        mx = (ins[i][1] > mx ? ins[i][1] : mx);
    }
    if (mx >= G->V) resizeGraph(G, mx + 1);
    resizeSSSP(T, G->V);

    int localAff = 0;
    int* newDist = malloc(T->size * sizeof(int));
    int* newParent = malloc(T->size * sizeof(int));
    if (!newDist || !newParent) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    memcpy(newDist, T->dist, T->size * sizeof(int));
    memcpy(newParent, T->parent, T->size * sizeof(int));

    #pragma omp parallel for schedule(static) reduction(+:localAff)
    for (int i = 0; i < insSz; i++) {
        int u = ins[i][0], v = ins[i][1], w = wgt[i];
        if (!edgeExists(G, u, v)) {
            #pragma omp critical
            addEdge(G, u, v, w);
        }
        int nd = (T->dist[v] != INF && T->dist[v] + w < newDist[u]) ? T->dist[v] + w : INF;
        if (nd < newDist[u]) {
            newDist[u] = nd;
            newParent[u] = v;
            T->aff[u] = 1;
            localAff++;
        }
        nd = (T->dist[u] != INF && T->dist[u] + w < newDist[v]) ? T->dist[u] + w : INF;
        if (nd < newDist[v]) {
            newDist[v] = nd;
            newParent[v] = u;
            T->aff[v] = 1;
            localAff++;
        }
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < T->size; i++) {
        T->dist[i] = newDist[i];
        T->parent[i] = newParent[i];
    }
    free(newDist); free(newParent);
    *affCount += localAff;
    return localAff > 0;
}

void propagateDeletion(const Graph* G, SSSP* T) {
    int more = 1;
    while (more) {
        more = 0;
        #pragma omp parallel for schedule(dynamic) reduction(||:more)
        for (int v = 0; v < G->V; v++) {
            if (!T->affDel[v]) continue;
            T->affDel[v] = 0;
            for (int i = 0; i < G->adj_sizes[v]; i++) {
                int wv = G->adj[v][i].to;
                if (T->parent[wv] == v) {
                    T->dist[wv] = INF;
                    T->parent[wv] = -1;
                    T->aff[wv] = 1;
                    T->affDel[wv] = 1;
                    more = 1;
                }
            }
        }
    }
}

int relaxLoop(const Graph* G, SSSP* T) {
    int changed = 0;
    int* newDist = malloc(T->size * sizeof(int));
    int* newParent = malloc(T->size * sizeof(int));
    if (!newDist || !newParent) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    memcpy(newDist, T->dist, T->size * sizeof(int));
    memcpy(newParent, T->parent, T->size * sizeof(int));

    #pragma omp parallel for schedule(dynamic) reduction(||:changed)
    for (int u = 0; u < G->V; u++) {
        if (!T->aff[u]) continue;
        T->aff[u] = 0;
        for (int i = 0; i < G->adj_sizes[u]; i++) {
            Edge e = G->adj[u][i];
            int nd = (T->dist[u] != INF) ? T->dist[u] + e.w : INF;
            if (nd < newDist[e.to]) {
                newDist[e.to] = nd;
                newParent[e.to] = u;
                T->aff[e.to] = 1;
                changed = 1;
            }
            nd = (T->dist[e.to] != INF) ? T->dist[e.to] + e.w : INF;
            if (nd < newDist[u]) {
                newDist[u] = nd;
                newParent[u] = e.to;
                T->aff[u] = 1;
                changed = 1;
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < T->size; i++) {
        T->dist[i] = newDist[i];
        T->parent[i] = newParent[i];
    }
    free(newDist); free(newParent);
    return changed;
}

/*------------ METIS CSR builder ------------*/
void buildCSR(const Graph* G, idx_t** xadj, idx_t** adjncy, idx_t** adjwgt, idx_t* nz) {
    *xadj = malloc((G->V + 1) * sizeof(idx_t));
    *nz = 0;
    for (int u = 0; u < G->V; u++) *nz += G->adj_sizes[u];
    *adjncy = malloc(*nz * sizeof(idx_t));
    *adjwgt = malloc(*nz * sizeof(idx_t));
    if (!*xadj || !*adjncy || !*adjwgt) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    idx_t ofs = 0;
    for (int u = 0; u < G->V; u++) {
        (*xadj)[u] = ofs;
        for (int i = 0; i < G->adj_sizes[u]; i++) {
            (*adjncy)[ofs] = G->adj[u][i].to;
            (*adjwgt)[ofs] = G->adj[u][i].w;
            ofs++;
        }
    }
    (*xadj)[G->V] = ofs;
}

double get_cpu_time() {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9; // seconds
}

double get_wall_time(struct timespec* ts) {
    struct timespec curr;
    clock_gettime(CLOCK_MONOTONIC, &curr);
    double sec = (curr.tv_sec - ts->tv_sec) + (curr.tv_nsec - ts->tv_nsec) / 1e9;
    *ts = curr;
    return sec;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Timing timing = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double start_cpu = get_cpu_time();
    double start_wall = MPI_Wtime();

    // Get hostname
    char hostname[256];
    gethostname(hostname, sizeof(hostname));

    if (argc < 4) {
        if (rank == 0) fprintf(stderr, "Usage: %s <edgeFile> <numEdits> <percentAdd>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    const char* edgeFile = argv[1];
    int numEdits = atoi(argv[2]);
    int pctAdd = atoi(argv[3]);
    if (pctAdd < 0 || pctAdd > 100) {
        if (rank == 0) fprintf(stderr, "percentAdd must be 0-100\n");
        MPI_Finalize();
        return 1;
    }
    int numIns = (int)(numEdits * pctAdd / 100.0 + 0.5);
    int numDel = numEdits - numIns;
    if (rank == 0) printf("Total edits: %d  Ins: %d  Del: %d\n", numEdits, numIns, numDel);

    Graph* G_global = initGraph(0);
    if (rank == 0) {
        FILE* fp = fopen(edgeFile, "r");
        if (!fp) { perror("fopen"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int u, v, w;
        while (fscanf(fp, "%d %d %d", &u, &v, &w) == 3) {
            if (w <= 0) {
                fprintf(stderr, "Negative or zero weight detected\n");
                fclose(fp);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            addEdge(G_global, u, v, w);
        }
        fclose(fp);
    }

    double comm_start = MPI_Wtime();
    int V = G_global->V;
    MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    if (rank != 0) resizeGraph(G_global, V);

    int E = 0;
    int (*triples)[3] = NULL;
    if (rank == 0) {
        for (int u = 0; u < V; u++)
            for (int i = 0; i < G_global->adj_sizes[u]; i++)
                if (u < G_global->adj[u][i].to) E++;
        triples = malloc(E * 3 * sizeof(int));
        if (!triples) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int idx = 0;
        for (int u = 0; u < V; u++)
            for (int i = 0; i < G_global->adj_sizes[u]; i++) {
                int v2 = G_global->adj[u][i].to;
                if (u < v2) {
                    triples[idx][0] = u;
                    triples[idx][1] = v2;
                    triples[idx][2] = G_global->adj[u][i].w;
                    idx++;
                }
            }
    }
    comm_start = MPI_Wtime();
    MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        triples = malloc(E * 3 * sizeof(int));
        if (!triples) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    }
    MPI_Bcast(triples, E * 3, MPI_INT, 0, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    if (rank != 0)
        for (int i = 0; i < E; i++)
            addEdge(G_global, triples[i][0], triples[i][1], triples[i][2]);
    if (rank == 0) printf("Loaded %d vertices, %d edges\n", V, E);

    idx_t* part = malloc(V * sizeof(idx_t));
    if (!part) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    if (size > 1) {
        idx_t *xadj, *adjncy, *adjwgt, nz;
        if (rank == 0) {
            buildCSR(G_global, &xadj, &adjncy, &adjwgt, &nz);
            idx_t nv = V, ncon = 1, objval;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            real_t* tpwgts = malloc(size * ncon * sizeof(real_t));
            if (!tpwgts) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
            for (int i = 0; i < size * ncon; i++) tpwgts[i] = 1.0f / size;
            real_t ubval = 1.05f, *ubvec = &ubval;
            int ret = METIS_PartGraphKway(&nv, &ncon, xadj, adjncy, NULL, NULL, adjwgt,
                                          &size, tpwgts, ubvec, options, &objval, part);
            if (ret != METIS_OK) {
                fprintf(stderr, "METIS partitioning failed: %d\n", ret);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            free(tpwgts); free(xadj); free(adjncy); free(adjwgt);
        }
        comm_start = MPI_Wtime();
        MPI_Bcast(part, V, MPI_INT, 0, MPI_COMM_WORLD);
        timing.comm_time += MPI_Wtime() - comm_start;
    } else {
        for (int i = 0; i < V; i++) part[i] = 0;
    }

    int* g2l = malloc(V * sizeof(int));
    int* l2g = NULL;
    if (!g2l) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    memset(g2l, -1, V * sizeof(int));
    int localV = 0;
    for (int i = 0; i < V; i++)
        if (part[i] == rank) localV++;
    l2g = malloc(localV * sizeof(int));
    if (!l2g) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    int local_idx = 0;
    for (int i = 0; i < V; i++)
        if (part[i] == rank) {
            g2l[i] = local_idx;
            l2g[local_idx] = i;
            local_idx++;
        }
    Graph* G = initGraph(localV);

    int* ghostG = NULL;
    int* ghostL = NULL;
    int ghostSz = 0;
    {
        char* seen = calloc(V, 1);
        if (!seen) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        for (int u = 0; u < V; u++)
            if (part[u] == rank)
                for (int i = 0; i < G_global->adj_sizes[u]; i++) {
                    int v2 = G_global->adj[u][i].to;
                    if (part[v2] != rank && !seen[v2]) { seen[v2] = 1; ghostSz++; }
                }
        ghostG = malloc(ghostSz * sizeof(int));
        ghostL = malloc(ghostSz * sizeof(int));
        if (!ghostG || !ghostL) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int idx = 0;
        for (int u = 0; u < V; u++)
            if (seen[u]) {
                ghostG[idx] = u;
                ghostL[idx] = localV + idx;
                g2l[u] = localV + idx;
                idx++;
            }
        free(seen);
    }
    resizeGraph(G, localV + ghostSz);
    for (int u = 0; u < V; u++)
        if (part[u] == rank)
            for (int i = 0; i < G_global->adj_sizes[u]; i++) {
                Edge e = G_global->adj[u][i];
                if (g2l[e.to] >= 0) addEdge(G, g2l[u], g2l[e.to], e.w);
            }

    SSSP* globalT = initSSSP(V);
    if (rank == 0) {
        dijkstra(G_global, globalT, 0);
    }
    comm_start = MPI_Wtime();
    MPI_Bcast(globalT->dist, V, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(globalT->parent, V, MPI_INT, 0, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    SSSP* localT = initSSSP(localV + ghostSz);
    for (int u = 0; u < V; u++)
        if (g2l[u] >= 0) {
            localT->dist[g2l[u]] = globalT->dist[u];
            localT->parent[g2l[u]] = globalT->parent[u] == -1 ? -1 : g2l[globalT->parent[u]];
        }

    int (*allDel)[2] = malloc(numDel * 2 * sizeof(int));
    int (*allIns)[2] = malloc(numIns * 2 * sizeof(int));
    int* allW = malloc(numIns * sizeof(int));
    if (!allDel || !allIns || !allW) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    if (rank == 0) {
        srand(12345);
        int di = 0, ii = 0;
        while (di < numDel) {
            int u = rand() % V, v = rand() % V;
            if (u != v && edgeExists(G_global, u, v)) {
                allDel[di][0] = u; allDel[di][1] = v; di++;
            }
        }
        while (ii < numIns) {
            int u = rand() % V, v = rand() % V;
            if (u != v && !edgeExists(G_global, u, v)) {
                allIns[ii][0] = u; allIns[ii][1] = v;
                allW[ii] = rand() % 10 + 1; ii++;
            }
        }
    }
    comm_start = MPI_Wtime();
    MPI_Bcast(allDel, numDel * 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(allIns, numIns * 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(allW, numIns, MPI_INT, 0, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;

    int delSz = 0, insSz = 0;
    int (*delLoc)[2] = malloc(numDel * 2 * sizeof(int));
    int (*insLoc)[2] = malloc(numIns * 2 * sizeof(int));
    int* insW = malloc(numIns * sizeof(int));
    if (!delLoc || !insLoc || !insW) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (int i = 0; i < numDel; i++) {
        int u = allDel[i][0], v = allDel[i][1];
        if (g2l[u] >= 0 && g2l[v] >= 0) {
            delLoc[delSz][0] = g2l[u];
            delLoc[delSz][1] = g2l[v];
            delSz++;
        }
    }
    for (int i = 0; i < numIns; i++) {
        int u = allIns[i][0], v = allIns[i][1];
        if (g2l[u] >= 0 && g2l[v] >= 0) {
            insLoc[insSz][0] = g2l[u];
            insLoc[insSz][1] = g2l[v];
            insW[insSz] = allW[i];
            insSz++;
        }
    }

    struct timespec t1, t2, t3, t4, t5;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    int affCount = 0;
    int anyChanges = 0;

    if (numDel > 0.5 * numEdits && rank == 0) {
        printf("Deletion-heavy workload, proceeding with updates\n");
    }

    for (int i = 0; i < delSz; i += BATCH_SIZE) {
        int currSz = (i + BATCH_SIZE <= delSz) ? BATCH_SIZE : (delSz - i);
        anyChanges |= processDeletions(G, localT, &delLoc[i], currSz, &affCount);
    }
    timing.del_time = get_wall_time(&t1);
    clock_gettime(CLOCK_MONOTONIC, &t2);

    
    for (int i = 0; i < insSz; i += BATCH_SIZE) {
        int currSz = (i + BATCH_SIZE <= insSz) ? BATCH_SIZE : (insSz - i);
        anyChanges |= processInsertions(G, localT, &insLoc[i], currSz, &insW[i], &affCount);
    }
    timing.ins_time = get_wall_time(&t2);
    clock_gettime(CLOCK_MONOTONIC, &t3);

    if (anyChanges) propagateDeletion(G, localT);
    timing.prop_time = get_wall_time(&t3);
    clock_gettime(CLOCK_MONOTONIC, &t4);

    int globalAff;
    comm_start = MPI_Wtime();
    MPI_Allreduce(&affCount, &globalAff, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    if (globalAff > 0.75 * V && rank == 0) {
        printf("Too many affected vertices (%d/%d), proceeding with updates\n", globalAff, V);
    }

    int* send_counts = malloc(size * sizeof(int));
    int* send_displs = malloc(size * sizeof(int));
    int* recv_counts = malloc(size * sizeof(int));
    int* recv_displs = malloc(size * sizeof(int));
    if (!send_counts || !send_displs || !recv_counts || !recv_displs) {
        fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1);
    }
    comm_start = MPI_Wtime();
    MPI_Allgather(&ghostSz, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    int total_ghosts = 0;
    recv_displs[0] = 0;
    for (int i = 0; i < size; i++) {
        send_counts[i] = (i == rank) ? ghostSz : 0;
        total_ghosts += recv_counts[i];
        if (i > 0) recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
    }
    int* all_ghostG = malloc(total_ghosts * sizeof(int));
    if (!all_ghostG) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    comm_start = MPI_Wtime();
    MPI_Allgatherv(ghostG, ghostSz, MPI_INT, all_ghostG, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;

    int globalCh = anyChanges;
    while (globalCh) {
        int localCh = relaxLoop(G, localT);
        int* send_dist = malloc(ghostSz * sizeof(int));
        int* send_parent = malloc(ghostSz * sizeof(int));
        int* recv_dist = malloc(total_ghosts * sizeof(int));
        int* recv_parent = malloc(total_ghosts * sizeof(int));
        if (!send_dist || !send_parent || !recv_dist || !recv_parent) {
            fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < ghostSz; i++) {
            send_dist[i] = localT->dist[ghostL[i]];
            send_parent[i] = localT->parent[ghostL[i]];
        }
        comm_start = MPI_Wtime();
        MPI_Allgatherv(send_dist, ghostSz, MPI_INT, recv_dist, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(send_parent, ghostSz, MPI_INT, recv_parent, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);
        timing.comm_time += MPI_Wtime() - comm_start;
        for (int i = 0; i < ghostSz; i++) {
            int global_v = ghostG[i];
            for (int j = 0; j < total_ghosts; j++) {
                if (all_ghostG[j] == global_v && recv_dist[j] < localT->dist[ghostL[i]]) {
                    localT->dist[ghostL[i]] = recv_dist[j];
                    localT->parent[ghostL[i]] = recv_parent[j];
                    localT->aff[ghostL[i]] = 1;
                    localCh = 1;
                }
            }
        }
        free(send_dist); free(send_parent); free(recv_dist); free(recv_parent);
        comm_start = MPI_Wtime();
        MPI_Allreduce(&localCh, &globalCh, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        timing.comm_time += MPI_Wtime() - comm_start;
    }
    timing.relax_time = get_wall_time(&t4);
    free(all_ghostG); free(send_counts); free(send_displs); free(recv_counts); free(recv_displs);

    /* 11) Output timings */
    timing.cpu_time = get_cpu_time() - start_cpu;
    comm_start = MPI_Wtime();
    MPI_Reduce(&timing.cpu_time, &timing.total_cpu_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    timing.comm_time += MPI_Wtime() - comm_start;
    timing.wall_time = MPI_Wtime() - start_wall;

    if (rank == 0) {
        printf("Phase timings (rank 0, host %s):\n", hostname);
        printf("  Deletion:   %.3f s\n", timing.del_time);
        printf("  Insertion:  %.3f s\n", timing.ins_time);
        printf("  Propagate:  %.3f s\n", timing.prop_time);
        printf("  Relax+Ghost:%.3f s\n", timing.relax_time);
        printf("Total communication time (rank 0): %.3f s\n", timing.comm_time);
        printf("Total CPU time (rank 0): %.3f s\n", timing.cpu_time);
        printf("Total wall clock time (rank 0): %.3f s\n", timing.wall_time);
        printf("Total cumulative CPU time (all ranks): %.3f s\n", timing.total_cpu_time);
    }
    for (int r = 0; r < size; r++) {
        if (r == rank) {
            printf("Rank %d (host %s):\n", rank, hostname);
            printf("  Comm time: %.3f s\n", timing.comm_time);
            printf("  CPU time:  %.3f s\n", timing.cpu_time);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* 12) Output SSSP tree to log.txt and compute accuracy */
    int* global_dist = NULL;
    int* global_parent = NULL;
    int* send_buf = malloc((localV + ghostSz) * 2 * sizeof(int));
    if (!send_buf) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (int i = 0; i < localV + ghostSz; i++) {
        send_buf[i*2] = localT->dist[i];
        send_buf[i*2 + 1] = localT->parent[i] == -1 ? -1 :
                            (i < localV ? l2g[localT->parent[i]] : ghostG[localT->parent[i] - localV]);
    }

    if (size == 1) {
        // Single-rank case: directly use localT data
        if (rank == 0) {
            global_dist = malloc(V * sizeof(int));
            global_parent = malloc(V * sizeof(int));
            if (!global_dist || !global_parent) {
                fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1);
            }
            for (int i = 0; i < V; i++) {
                global_dist[i] = localT->dist[i];
                global_parent[i] = localT->parent[i];
            }
        }
    } else {
        // Multi-rank case: use MPI_Gatherv to collect data
        int* recv_counts = malloc(size * sizeof(int));
        int* recv_displs = malloc(size * sizeof(int));
        if (!recv_counts || !recv_displs) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int send_count = (localV + ghostSz) * 2;
        comm_start = MPI_Wtime();
        MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
        timing.comm_time += MPI_Wtime() - comm_start;

        int total_recv = 0;
        recv_displs[0] = 0;
        for (int i = 0; i < size; i++) {
            total_recv += recv_counts[i];
            if (i > 0) recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        }

        int* recv_buf = NULL;
        if (rank == 0) {
            recv_buf = malloc(total_recv * sizeof(int));
            global_dist = malloc(V * sizeof(int));
            global_parent = malloc(V * sizeof(int));
            if (!recv_buf || !global_dist || !global_parent) {
                fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1);
            }
            for (int i = 0; i < V; i++) {
                global_dist[i] = INF;
                global_parent[i] = -1;
            }
        }

        comm_start = MPI_Wtime();
        MPI_Gatherv(send_buf, send_count, MPI_INT, recv_buf, recv_counts, recv_displs, MPI_INT, 0, MPI_COMM_WORLD);
        timing.comm_time += MPI_Wtime() - comm_start;

        if (rank == 0) {
            for (int r = 0; r < size; r++) {
                int r_localV = 0;
                for (int i = 0; i < V; i++)
                    if (part[i] == r) r_localV++;
                int r_ghostSz = recv_counts[r] / 2 - r_localV;
                int offset = recv_displs[r];
                for (int i = 0; i < V; i++) {
                    if (part[i] == r) {
                        int li = g2l[i];
                        global_dist[i] = recv_buf[offset + li*2];
                        global_parent[i] = recv_buf[offset + li*2 + 1];
                    }
                }
            }
            free(recv_buf);
        }
        free(recv_counts);
        free(recv_displs);
    }
    free(send_buf);

    /* Compute accuracy and log SSSP tree */
    double accuracy = 0.0;
    if (rank == 0) {
        // Apply edits to G_global for ground-truth
        for (int i = 0; i < numDel; i++) {
            int u = allDel[i][0], v = allDel[i][1];
            if (edgeExists(G_global, u, v)) removeEdge(G_global, u, v);
        }
        for (int i = 0; i < numIns; i++) {
            int u = allIns[i][0], v = allIns[i][1], w = allW[i];
            if (!edgeExists(G_global, u, v)) addEdge(G_global, u, v, w);
        }

        // Compute source's connected component
        char* reachable = malloc(V * sizeof(char));
        if (!reachable) { fprintf(stderr, "Memory allocation failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int size_reachable = bfs_component(G_global, 0, reachable);

        // Compute reference SSSP
        SSSP* refT = initSSSP(V);
        dijkstra(G_global, refT, 0);

        // Compute accuracy
        int correct_reachable = 0, correct_unreachable = 0;
        for (int i = 0; i < V; i++) {
            if (reachable[i]) {
                if (global_dist[i] == refT->dist[i]) correct_reachable++;
            } else {
                if (global_dist[i] == INF) correct_unreachable++;
            }
        }
        accuracy = (V > 0) ? ((correct_reachable + correct_unreachable) * 100.0 / V) : 0.0;
        free(reachable);
        freeSSSP(refT);

        // Log SSSP tree and accuracy to log.txt
        FILE* fp = fopen("log.txt", "w");
        if (!fp) { perror("fopen log.txt"); MPI_Abort(MPI_COMM_WORLD, 1); }
        fprintf(fp, "Final SSSP Tree:\nVertex\tDistance\tParent\n");
        for (int i = 0; i < V; i++) {
            fprintf(fp, "%d\t%d\t\t%d\n", i, global_dist[i], global_parent[i]);
        }
        fclose(fp);
    }

    /* 13) Cleanup */
    freeGraph(G_global);
    if (triples) free(triples);
    freeGraph(G);
    freeSSSP(globalT);
    freeSSSP(localT);
    free(part);
    free(g2l);
    free(l2g);
    free(delLoc);
    free(insLoc);
    free(insW);
    free(allDel);
    free(allIns);
    free(allW);
    free(ghostG);
    free(ghostL);
    if (rank == 0) {
        free(global_dist);
        free(global_parent);
    }

    MPI_Finalize();
    return 0;
}