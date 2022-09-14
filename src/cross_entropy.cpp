#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cross_entropy.h"

static const double e = 2.71828;


int main(int argc, char* argv[]){


    //Graph Variables
    char* filename;
    int N;
    int max_degree;

    int* degrees;
    int** neighbours;
    int** edges;

    //Parameters
    int n = 50;
    int m = 10;
    int r = 5;
    double rho = 0.1;
    double alpha = 0.5;


    //Calculated vars
    double delta;

    int* dombest;
    int** domsets;
    double* P;
    double* Pstar;
    double* L;

    int label_offset = 1;
    int seed;
    int output_types = 0;
    
    bool has_file = 0;

    bool (*dom_func)(int*, int, int*, int**) = &dominates;

    for (int i = 0; i < argc; i++)
    {
        char* s = argv[i];
        
        if(!strcmp(s, "-f")){
            filename = argv[++i];
            has_file = 1;

            if(i + 1 < argc && argv[i+1][0] != '-'){
                label_offset = atoi(argv[++i]);
            }
        } else if(!strcmp(s, "-n")){
            n = atoi(argv[++i]);
        } else if(!strcmp(s, "-m")){
            m = atoi(argv[++i]);
        } else if(!strcmp(s, "-r")){
            r = atoi(argv[++i]);
        } else if(!strcmp(s, "-R")){
            rho = atof(argv[++i]);
        } else if(!strcmp(s, "-a")){
            alpha = atof(argv[++i]);
        } else if(!strcmp(s, "-s")){
            seed = atof(argv[++i]);
        } else if(!strcmp(s, "-o")){
            output_types = atof(argv[++i]);
        } else if(!strcmp(s, "-d")){
            s = argv[++i];
            if(!strcmp(s, "d")){
                dom_func = &dominates;
            } else if(!strcmp(s, "s")){
                dom_func = &secure_dominates;
            }

        }
    }
    
    if(!has_file){
        fprintf(stdout, "Please provide a edge list file with -f \"filename\"");
        return 1;
    }

    srand(seed);

    int edge_count = read_edges(edges, filename, label_offset);


    make_graph(N, degrees, neighbours, edge_count, edges);
    
    

    dombest = (int*) malloc(N*sizeof(int));//new int[N];
    P = (double*) malloc(N*sizeof(double));//new double[N];
    Pstar = (double*) malloc(N*sizeof(double));//new double[N];
    for (int i = 0;i < N;i++)
    {
        dombest[i] = 1;
        P[i] = 1.0/double(N);
        Pstar[i] = 0;
    }

    domsets = (int**) malloc(n*sizeof(int*));//new int*[n];
    L =  (double*) malloc(n*sizeof(double));//new double[n];
    int t = 0;
    while(true){

        for (int i = 0; i < n; i++)
        {
            make_domset(domsets[i], N, degrees, neighbours, P, (*dom_func));
            L[i] = calculate_score(N, domsets[i]);
        }


        sort_domsets(L, domsets, n, N);
        
        if(calculate_score(N, dombest) > L[0]){
            dombest = domsets[0];
        }

        if(t > r){
            break;
        }

        delta = -1*L[0]/log(rho);

        double psum = 0;
        for (int i = 0; i < N; i++)
        {
            Pstar[i] = calculate_Pstar(i, m, L, domsets, delta);
            psum += Pstar[i];
        }

        for (int i = 0; i < N; i++)
        {
            P[i] = (1-alpha)*P[i] + alpha*(Pstar[i]/psum);
        }

        t++;
    
    }

    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += dombest[i];
    }
            

    fprintf(stdout, "Best is %d guards\n", sum);

    for (int i  = 0; i < N; i++)
    {
        if(output_types == 1){
            fprintf(stdout, "%d ", dombest[i]);
        } else if(output_types == 2 && dombest[i]){
            fprintf(stdout, "%d ", i+label_offset);
        } else if(output_types == 3){
            fprintf(stdout, "%d    ", dombest[i]);
        }
        
    }
    if(output_types > 0){
        fprintf(stdout, "\n");
    }
    if(output_types == 3){
        for (int i  = 0; i < N; i++)
        {
        
            fprintf(stdout, "%.2f ", P[i]);
        
        
        }
    
        fprintf(stdout, "\n");
    }




    for (int i = 0; i < edge_count; i++)
    {
        free(edges[i]);
    }
    for(int i = 0; i < N; i++)
    {
        free(neighbours[i]);
    }
    for (int i = 0; i < n; i++)
    {
        free(domsets[i]);
    }

    free(dombest);
    free(P);
    free(Pstar);
    free(L);
    free(degrees);
    free(edges);
    free(neighbours);
    free(domsets);
    
    

    return 0;
}

/**
 * @brief Reads a list of edges from a given file. 
 * 
 * @param edges will contain list of edges
 * @param filename name of file to open
 * @param label_offset either 0 or 1 depending on indexing of vertices
 * @return -1 if no file, number of edges otherwise.
 */
int read_edges(int**& edges, char* filename, int label_offset){
    //std::ifstream file(filename);
    FILE* file = fopen(filename, "r");

    if(!file){
        fprintf(stdout, "File does not exist\n");
        return -1;
    }


    int** edge_buffer = (int**) malloc(10000*sizeof(int*));//new int*[10000];

    int count = 0;

    while(true){
        int a = 0;
        int b = 0;
        int check = fscanf(file, " %d %d ", &a, &b);
        if(check == EOF) break;

        
        edge_buffer[count] = (int*) malloc(2*sizeof(int));//new int[2];
        edge_buffer[count][0] = a - label_offset;
        edge_buffer[count++][1] = b - label_offset;
        
    }

    edges = (int**) malloc(count*sizeof(int*));//new int*[count];
    for (int i = 0; i < count; i++)
    {
        edges[i] = (int*) malloc(2*sizeof(int));//new int[2];
        edges[i][0] = edge_buffer[i][0];
        edges[i][1] = edge_buffer[i][1];
        free(edge_buffer[i]);
    }
    
    free(edge_buffer);

    
    return count;

}

/**
 * @brief Creates a list of neighbours given a set of edges
 * 
 * @param N will be number of vertices in the graph
 * @param degrees will belist of degrees for each vertex
 * @param neighbours will be list of neighbours for each vertex
 * @param edge_count number of edges
 * @param edges list of edges
 */
void make_graph(int &N, int* &degrees, int** &neighbours, int edge_count, int** edges){
    N = -1;
    for (int i = 0; i < edge_count; i++)
    {
        if(edges[i][0] > N){
            N = edges[i][0];
        }
        if(edges[i][1] > N){
            N = edges[i][1];
        }
    }
    N++;

    degrees = (int*) malloc(N*sizeof(int));//new int[N];
    memset(degrees, 0, N*sizeof(int));
    
    for (int i = 0; i < edge_count; i++)
    {
        int a = edges[i][0];
        int b = edges[i][1];

        degrees[a]++;
        degrees[b]++;
    }

    

    neighbours = (int**) malloc(N*sizeof(int*));//new int*[N];
    for (int i = 0; i < N; i++)
    {
        neighbours[i] = (int*) malloc(degrees[i]*sizeof(int));//new int[degrees[i]];
    }
    
    int counter[N] = {0};
    for (int i = 0; i < edge_count; i++)
    {

        int a = edges[i][0];
        int b = edges[i][1];


        neighbours[a][counter[a]++] = b;
        neighbours[b][counter[b]++] = a;

    }
    
}

/**
 * @brief Fills a given array with a dominating set
 * 
 * @param domset array that will contain a domset
 * @param N number of vertices for the domset
 * @param degrees degree of each vertex in the graph to be dominated
 * @param neighbours neighbours of each vertex in the graph to be dominated
 * @param P probabilities that each vertex will be selected for the domset
 */
void make_domset(int* &domset, int N, int* degrees, int** neighbours, double* P, bool (*dom_func)(int*, int, int*, int**)){
    domset = (int*) malloc(N*sizeof(int));//new int[N];
    memset(domset, 0, N*sizeof(int));

    double Ptemp[N];
    memcpy(Ptemp, P, sizeof(*P)*N);

    while(!(*dom_func)(domset, N, degrees, neighbours)){
        int ind = weight_rand(N, Ptemp);
        domset[ind] = 1;
        Ptemp[ind] = 0;
    }

    for (int i = 0; i < N; i++){
        if(domset[i]){
            domset[i] = 0;
            if(!(*dom_func)(domset, N, degrees, neighbours)){
                domset[i] = 1;
            }
        }
    }


}

/**
 * @brief Calculates the score of a given domset
 * 
 * @param N number of vertices in domset
 * @param domset the domset
 * @return the score of the domset, lower is better
 */
double calculate_score(int N, int* domset){
    int sum = 0;

    for (int i = 0; i < N; i++)
    {
        sum += domset[i];
    }
    
    return sum;
}

/**
 * @brief Sorts both L and domsets based on the values in L
 * 
 * @param L unsorted list of scores. will be sorted
 * @param domsets unsorted list of dominating sets. will be sorted
 * @param n number of domsets
 * @param N number of vertices for the domsets
 */
void sort_domsets(double* &L, int** &domsets, int n, int N){
    double to_sort[n][2];
    for (int i = 0; i < n; i++)
    {
        to_sort[i][0] = i;
        to_sort[i][1] = L[i];
    }
    
    qsort(to_sort, n, sizeof(*to_sort), compare_scores);
    
    double Ltemp[n];
    int domset_temp[n][N];

    for (int i = 0; i < n; i++)
    {
        Ltemp[i] = L[(int)to_sort[i][0]];
        for (int j = 0; j < N; j++)
        {
            domset_temp[i][j] = domsets[(int)to_sort[i][0]][j];
        }
        
    }

    for (int i = 0; i < n; i++)
    {
        L[i] = Ltemp[i];
        for (int j = 0; j < N; j++)
        {
            domsets[i][j] = domset_temp[i][j];
        }
        
    }
    
}

/**
 * @brief Compares the elements of a list of pairs, based on the second value.
 * 
 * @param a first pair
 * @param b second pair
 * @return 0 if they are equal, -1 if a < b and 1 if a > b
 */
int compare_scores( const void* a, const void* b)
{
    double L_a = ((double*)a)[1];
    double L_b = ((double*)b)[1];

     if ( L_a == L_b ) return 0;
     else if ( L_a < L_b ) return -1;
     else return 1;
}

/**
 * @brief Calculates the updated probability distribution for a given index.
 * 
 * @param i index of the vertex
 * @param m number of best domsets to choose from
 * @param L list of scores for each domset, sorted low to high
 * @param domsets list of domsets, sorted by score
 * @param delta pre-calculated value for scoring
 * @return the updated probability for the given index.
 */
double calculate_Pstar(int i, int m, double* L, int** domsets, double delta){

    double sum_num = 0;
    double sum_den = 0;
    
    for (int j = 0; j < m; j++)
    {
        double val = pow(e, -L[j]/delta);
        sum_den += val;
        if(domsets[j][i]){
            sum_num += val;
        }
    }
    
    return sum_num/sum_den;
}

/**
 * @brief Determine if the given dominating set dominates the graph defined by the list of neighbours.
 * 
 * @param domset The dominating set to check
 * @param N Number of vertices in the graph
 * @param degrees list containing the degree of each vertex
 * @param neighbours the neighbours of each vertex
 * @return true if the given set is dominating
 * @return false if the given set is not dominating
 */
bool dominates(int* domset, int N, int* degrees, int** neighbours){
    int dommed[N] = {0};

    for (int i = 0; i < N; i++)
    {
        if(domset[i]){
            dommed[i] = 1;
            for (int j = 0; j < degrees[i]; j++)
            {
                dommed[neighbours[i][j]] = 1;
            }
            
        }
    }

    for (int i = 0; i < N; i++)
    {
        if(!dommed[i]){
            return false;
        }
    }
    
    
    return true;
}

bool secure_dominates(int* domset, int N, int* degrees, int** neighbours){
    if(!dominates(domset, N, degrees, neighbours)){
        return false;
    }

    for (int i = 0; i < N; i++)
    {
        if (domset[i] == 0)
        {
            bool dommed = false;
            for (int j = 0; j < degrees[i]; j++)
            {
                int nb = neighbours[i][j];

                if(domset[nb]){
                    domset[i]++;
                    domset[nb]--;
                    dommed = dominates(domset, N, degrees, neighbours);
                    domset[i]--;
                    domset[nb]++;

                    if(dommed){
                        continue;
                    }
                }
            }
            if(!dommed){
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Given a distribution P of length N, randomly select an index based on the distribution.
 * 
 * @param N number of option in P
 * @param P some distribution of values
 * @return index of the selected value
 */
int weight_rand(int N, double* P){

    double sum = P[0];
    double cumulative[N] = {sum};
    for (int i = 1; i < N; i++)
    {
        sum += P[i];
        cumulative[i] = P[i] + cumulative[i-1];
        
    }

    double val = sum*rand()/double(RAND_MAX);

    for (int i = 0; i < N; i++)
    {
        if(val < cumulative[i]){
            return i;
        }
    }

    return 0;
}