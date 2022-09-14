int read_edges(int**&, char*, int);
void make_graph(int&, int*&, int**&, int, int**);
void make_domset(int*&, int, int*, int**, double*, bool (*)(int*,int,int*,int**));
double calculate_score(int, int*);
void sort_domsets(double*&, int**&, int, int);
int compare_scores(const void*, const void*);
double calculate_Pstar(int, int, double*, int**, double);
bool dominates(int*, int, int*, int**);
bool total_dominates(int*, int, int*, int**);
bool two_dominates(int*, int, int*, int**);
bool secure_dominates(int*, int, int*, int**);
int weight_rand(int, double*);


