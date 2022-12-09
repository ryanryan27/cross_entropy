typedef bool (*domfunc)(struct DomUpdater&, int, int*, struct Graph);

enum variant{
    domination,
    total,
    two,
    secure,
    connected
};

struct Params {
    char* filename;
    char* out_file;
    int n = 50;
    int m = 10;
    int r = 5;
    double rho = 0.1;
    double alpha = 0.5;
    char* dom_type_str = (char*)"d";
    variant dom_type = domination;

    int label_offset = 1;
    int seed = 0;
    int iterations = 1;
    int output_types = 2;
    double timeout = -1;

};

struct Graph {
    int N;
    int M;
    int max_degree;
    int* degrees;
    int** neighbours;
    int** edges;

    int* degrees3;
    int** three_aparts;
};

struct CEUpdater {
    double delta;
    int* best_domset_this_iteration;
    int** domsets;
    double* P;
    double* Pstar;
    double* L;

    int* results;
    int best;
    int* best_domset_overall;
    int best_seed;
    int domset_possible;
    bool timed_out;
    int loops_without_change;
    clock_t start;
    double total_time;
};

struct DomUpdater {
    int* dommed;
    int domsum;

    int* secure_dommed;
    int** secure_dom_neighbours;

    bool subtracting;
    int num_components;
    int* components;

    double* Ptemp;
    double sumP;

    int* links;

    domfunc dom_func;
};

void print_instructions();
void run_cross_entropy(Params);
void cross_entropy_main_loop(CEUpdater&, Graph, Params);
void handle_params(Params&, int, char*[]);
int read_edges(Graph&, Params);
void make_graph(Graph&, Params);
void init_updater(CEUpdater&, Graph, Params);
void init_dom_updater(DomUpdater&, Graph, Params, CEUpdater);
int make_domset(int*&, Graph, CEUpdater, Params);
double calculate_score(int, int*);
void sort_domsets(CEUpdater&, Graph, Params);
int compare_scores(const void*, const void*);
double calculate_Pstar(int, CEUpdater, Params);
void set_domfunc(DomUpdater&, Params);
void print_output(CEUpdater, Params, Graph);
void destruct_memory(Graph, CEUpdater, Params);
void destruct_memory(DomUpdater, Graph, Params);

bool dominates(DomUpdater&, int, int*, Graph);
bool total_dominates(DomUpdater&, int, int*, Graph);
bool two_dominates(DomUpdater&, int, int*, Graph);
bool secure_dominates(DomUpdater&, int, int*, Graph);
bool connected_dominates(DomUpdater&, int, int*, Graph);
bool is_connected(int*, Graph);
bool is_connected(DomUpdater&, int, int*, Graph);

int weight_rand_acc(int, double*, double, int*&);

void generate_three_apart(Graph& graph);
bool can_secure_dom(int, int, int*, Graph, DomUpdater);