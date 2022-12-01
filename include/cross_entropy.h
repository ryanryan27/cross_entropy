typedef bool (*domfunc)(struct DomUpdater&, int, int*, struct Graph);

struct Params {
    char* filename;
    int n = 50;
    int m = 10;
    int r = 5;
    double rho = 0.1;
    double alpha = 0.5;
    char* dom_type = (char*)"d";

    int label_offset = 1;
    int seed = 0;
    int iterations = 1;
    int output_types = 0;

};

struct Graph {
    int N;
    int M;
    int max_degree;
    int* degrees;
    int** neighbours;
    int** edges;
};

struct CEUpdater {
    double delta;
    int* dombest;
    int** domsets;
    double* P;
    double* Pstar;
    double* L;

    int* results;
    int best;
    int* best_domset;
};

struct DomUpdater {
    int* dommed;
    int domsum;

    double* Ptemp;
    double sumP;

    int* links;

    domfunc dom_func;
};


int read_edges(Graph&, Params);
void make_graph(Graph&);
void init_updater(CEUpdater&, Graph, Params);
int make_domset(int*&, Graph, CEUpdater, Params);
double calculate_score(int, int*);
void sort_domsets(CEUpdater&, Graph, Params);
int compare_scores(const void*, const void*);
double calculate_Pstar(int, CEUpdater, Params);
void set_domfunc(DomUpdater&, Params);

bool dominates(DomUpdater&, int, int*, Graph);
bool total_dominates(int*, int*&, int&, int, int, int*, int**);
bool two_dominates(int*, int*&, int&, int, int, int*, int**);
bool secure_dominates(int*, int, int*, int**);
bool connected_dominates(int*, int*&, int&, int, int, int*, int**);
bool connected(int*, int, int*, int**);

int weight_rand(int, double*, double);
int weight_rand_acc(int, double*, double, int*&);