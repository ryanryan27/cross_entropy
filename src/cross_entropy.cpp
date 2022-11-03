#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
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
    int seed = 0;
    int iterations = 1;
    int output_types = 0;

    double timeout = -1;
    bool out_of_time = false;
    
    bool has_file = 0;

    char* dom_type =(char*)"d";
    bool (*dom_func)(int*, int*&, int&, int, int, int*, int**) = &dominates;

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
        } else if(!strcmp(s, "-t")){
            timeout = atof(argv[++i]);
        } else if(!strcmp(s, "-o")){
            output_types = atof(argv[++i]);
        } else if(!strcmp(s, "-i")){
            iterations = atoi(argv[++i]);
        } else if(!strcmp(s, "-d")){
            dom_type = argv[++i];
            if(!strcmp(dom_type, "d")){
                dom_func = &dominates;
            } else if(!strcmp(dom_type, "s")){
                //dom_func = &secure_dominates;
            } else if(!strcmp(dom_type, "t")){
                dom_func = &total_dominates;
            } else if(!strcmp(dom_type, "2")){
                dom_func = &two_dominates;
            } else if(!strcmp(dom_type, "c")){
                dom_func = &connected_dominates;
            }

        }
    }
    
    if(!has_file){
        fprintf(stdout, "Please provide a edge list file with -f \"filename\"");
        return 1;
    }

    
    clock_t start = clock();

    int edge_count = read_edges(edges, filename, label_offset);


    make_graph(N, degrees, neighbours, edge_count, edges);

    

    dombest = (int*) malloc(N*sizeof(int));
    P = (double*) malloc(N*sizeof(double));
    Pstar = (double*) malloc(N*sizeof(double));


    domsets = (int**) malloc(n*sizeof(int*));
    L =  (double*) malloc(n*sizeof(double));

    for (int i = 0; i < n; i++)
    {
        domsets[i] = (int*) malloc(N*sizeof(int));
    }
    

    int results[iterations] = {0};
    int best = N;
    int best_domset[N] = {1};

    bool timed_out = false;
    int its_completed = 0;

    int domset_possible = 1;

    
    for (int i = 0; i < iterations; i++)
    {
        if (!domset_possible) break;
        
        srand(seed + i);

        for (int j = 0;j < N;j++)
        {
            dombest[j] = 1;
            P[j] = 1.0/double(N);
            Pstar[j] = 0;
        }

        int t = 0;

        while(true){
            
        
            for (int j = 0; j < n; j++)
            {
                domset_possible = make_domset(domsets[j], N, degrees, neighbours, P, (*dom_func));
            
                if (!domset_possible) break;

                L[j] = calculate_score(N, domsets[j]);
            }

        

            if (!domset_possible) break;

            sort_domsets(L, domsets, n, N);

            
            if(calculate_score(N, dombest) > L[0]){
                memcpy(dombest, domsets[0], N*sizeof(int));
                t = 0;
            }

            if(t > r){
                break;
            }

            if(timeout > 0 && ((double)(clock()-start)/CLOCKS_PER_SEC) > timeout){
                timed_out = true;
                break;
            }


            delta = -1*L[0]/log(rho);

            double psum = 0;
            for (int j = 0; j < N; j++)
            {
                Pstar[j] = calculate_Pstar(j, m, L, domsets, delta);
                psum += Pstar[j];
            }

            for (int j = 0; j < N; j++)
            {
                P[j] = (1-alpha)*P[j] + alpha*(Pstar[j]/psum);
            }

            
            t++;
        
        }

        int sum = 0;
        for (int j = 0; j < N; j++)
        {
            sum += dombest[j];
        }

        results[i] = sum;

        if(results[i] <best){
            best = results[i];
            for (int j = 0; j < N; j++)
            {
                best_domset[j] = dombest[j];
            }
            
        }

        its_completed = i+1;
        if(timed_out){
            
            break;
        }
    }

    double total_time  = (double)(clock() - start)/CLOCKS_PER_SEC;
    if(domset_possible){
        if(output_types > 0){
            fprintf(stdout, "Best is %d guards\n", best);
            fprintf(stdout, "Time taken: %0.3f \n", total_time);
            fprintf(stdout, "Dominating set: \n");
        }

        for (int i  = 0; i < N; i++)
        {
            if(output_types == 1){
                fprintf(stdout, "%d ", best_domset[i]);
            } else if(output_types == 2 && best_domset[i]){
                fprintf(stdout, "%d ", i+label_offset);
            } else if(output_types == 3){
                fprintf(stdout, "%d    ", best_domset[i]);
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

        if(output_types == -1){
            fprintf(stdout, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", filename, dom_type, n, m, r, rho, alpha, best, total_time);
        }
    } else {
        if(output_types == -1){
            fprintf(stdout, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", filename, dom_type, n, m, r, rho, alpha, -1, 0.00);
        } else {
            fprintf(stdout, "Unable to dominate graph - probably disconnected\n");
        }
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
    FILE* file = fopen(filename, "r");

    if(!file){
        fprintf(stdout, "File does not exist\n");
        return -1;
    }


    int** edge_buffer = (int**) malloc(10000*sizeof(int*));

    int count = 0;

    while(true){
        int a = 0;
        int b = 0;
        int check = fscanf(file, " %d %d ", &a, &b);
        if(check == EOF) break;

        
        edge_buffer[count] = (int*) malloc(2*sizeof(int));
        edge_buffer[count][0] = a - label_offset;
        edge_buffer[count++][1] = b - label_offset;
        
    }

    edges = (int**) malloc(count*sizeof(int*));
    for (int i = 0; i < count; i++)
    {
        edges[i] = (int*) malloc(2*sizeof(int));
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

    degrees = (int*) malloc(N*sizeof(int));
    memset(degrees, 0, N*sizeof(int));
    
    for (int i = 0; i < edge_count; i++)
    {
        int a = edges[i][0];
        int b = edges[i][1];

        degrees[a]++;
        degrees[b]++;
    }

    

    neighbours = (int**) malloc(N*sizeof(int*));
    for (int i = 0; i < N; i++)
    {
        neighbours[i] = (int*) malloc(degrees[i]*sizeof(int));
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
 * @return 1 if domset is possible, 0 if not
 */
int make_domset(int* &domset, int N, int* degrees, int** neighbours, double* P, bool (*dom_func)(int*, int*&, int&, int, int, int*, int**)){
    
    memset(domset, 0, N*sizeof(int));

    int* dommed = (int*) malloc(N*sizeof(int));
    memset(dommed, 0, N*sizeof(int));

    int domsum = 0;

    double Ptemp[N];
    memcpy(Ptemp, P, sizeof(*P)*N);

    double sumP = 0;
    

    for (int i = 0; i < N; i++)
    {
        Ptemp[i] += 1;
        sumP += Ptemp[i];
    }

    int* links = (int*) malloc((2*N+1)*sizeof(int));
    links[N-1] = -1;
    links[N] = -1;
    links[2*N] = 0;
    
    for(int i = 0; i < N-1; i++){
        links[i] = i+1;
        links[N + i + 1] = i;
    }
    
    // fprintf(stdout, "adding \n");

    int domcount = 0;
    int ind;
    do {
        if(domcount == N){
            free(dommed);
            //free(links);
            return 0;    
        }

        // ind = weight_rand(N, Ptemp, sumP);
        ind = weight_rand_acc(N, Ptemp, sumP, links);

        domset[ind] = 1;
        domcount++;
        sumP -= Ptemp[ind];
        Ptemp[ind] = 0;

    } while(!(*dom_func)(domset, dommed, domsum, ind, N, degrees, neighbours));


    for (int i = 0; i < N; i++){

        Ptemp[i] = 1.0 - P[i];

        if(Ptemp[i] < 0 ){
            fprintf(stdout, "oh no \n");
            exit(0);
        }

    }


    sumP = 0;

    for (int i = 0; i < N; i++)
    {
        sumP += Ptemp[i];
    }

    links[N-1] = -1;
    links[N] = -1;
    links[2*N] = 0;
    
    for(int i = 0; i < N-1; i++){
        links[i] = i+1;
        links[N + i + 1] = i;
    }

    //  fprintf(stdout, "removing \n");

    for (int i = 0; i < N; i++)
    {
        // int ind = weight_rand(N, Ptemp, sumP);
        int ind = weight_rand_acc(N, Ptemp, sumP, links);

        sumP -= Ptemp[ind];
        Ptemp[ind] = 0;

        if(domset[ind]){
            domset[ind] = 0;
            if(!(*dom_func)(domset, dommed, domsum, ind, N, degrees, neighbours)){
                domset[ind] = 1;
                (*dom_func)(domset, dommed, domsum, ind, N, degrees, neighbours);
            }
        }

    }

    //fprintf(stdout, "here \n");

    free(dommed);
    //free(links);
    return 1;

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
bool dominates(int* domset, int* &dommed, int &domsum, int added, int N, int* degrees, int** neighbours){

    if(domset[added]){
        if(!dommed[added]) domsum++;
        dommed[added] = 1;
        for (int j = 0; j < degrees[added]; j++)
        {
            if(!dommed[neighbours[added][j]]) 
            {
                domsum++;
                dommed[neighbours[added][j]] = 1;
            }
        }

    } else {

        bool undommed = true;
        for (int j = 0; j < degrees[added]; j++)
        {
            int nbr = neighbours[added][j];
            if(domset[nbr]){
                undommed = false;
                break;
            }
        }

        if(undommed){
            dommed[added] = 0;
            domsum--;
            //return false;
        } 

        for (int i = 0; i < degrees[added]; i++)
        {
            int neighbour = neighbours[added][i];

            if(domset[neighbour]) continue;

            undommed = true;
            for (int j = 0; j < degrees[neighbour]; j++)
            {
                int nbr = neighbours[neighbour][j];
                if(domset[nbr]){
                    undommed = false;
                    break;
                }
            }

            if(undommed){
                dommed[neighbour] = 0;
                domsum--;
                //return false;
            } 
        }
    }

    return domsum == N;

}

/**
 * @brief Determines if the given domset is total dominating. 
 * Maintains a list of domination status of each vertex.
 * 
 * @param domset 
 * @param dommed 
 * @param domsum 
 * @param added 
 * @param N 
 * @param degrees 
 * @param neighbours 
 * @return true 
 * @return false 
 */
bool total_dominates(int* domset, int* &dommed, int &domsum, int added, int N, int* degrees, int** neighbours){

    if(domset[added]){

        for (int j = 0; j < degrees[added]; j++)
        {
            if(!dommed[neighbours[added][j]]) 
            {
                domsum++;
                dommed[neighbours[added][j]] = 1;
            }
        }

    } else {

        for (int i = 0; i < degrees[added]; i++)
        {
            int neighbour = neighbours[added][i];

            bool undommed = true;
            for (int j = 0; j < degrees[neighbour]; j++)
            {
                int nbr = neighbours[neighbour][j];
                if(domset[nbr]){
                    undommed = false;
                    break;
                }
                
            }

            if(undommed){
                dommed[neighbour] = 0;
                domsum--;
            } 
        }
    }

    return domsum == N;

}


/**
 * @brief Determine if the given domset is two-dominating
 * 
 * @param domset 
 * @param dommed 
 * @param domsum 
 * @param added 
 * @param N 
 * @param degrees 
 * @param neighbours 
 * @return true 
 * @return false 
 */
bool two_dominates(int* domset, int* &dommed, int &domsum, int added, int N, int* degrees, int** neighbours){

    if(domset[added]){
        if(dommed[added] < 2) domsum++;
        dommed[added] = 2;
        for (int j = 0; j < degrees[added]; j++)
        {
            if(!domset[neighbours[added][j]]) 
            {
                dommed[neighbours[added][j]]++;

                if(dommed[neighbours[added][j]] == 2) domsum++;
                
            }
        }

    } else {

        int domset_nbrs = 0;
        for (int j = 0; j < degrees[added]; j++)
        {
            int nbr = neighbours[added][j];
            if(domset[nbr]){
                domset_nbrs++;
            }
        }
        dommed[added] = domset_nbrs;
        if(domset_nbrs < 2){
            domsum--;
            return false;
        } 


        for (int i = 0; i < degrees[added]; i++)
        {
            int neighbour = neighbours[added][i];

            if(domset[neighbour]) continue;

            domset_nbrs = 0;
            for (int j = 0; j < degrees[neighbour]; j++)
            {
                int nbr = neighbours[neighbour][j];
                if(domset[nbr]){
                    domset_nbrs++;
                }
            }
            
            dommed[neighbour] = domset_nbrs;
             
            if(domset_nbrs < 2){
                domsum--;
                return false;
            } 
        }
    }

    return domsum == N;

}


/**
 * @brief Determine if the given dominating set securely dominates the graph defined by the list of neighbours.
 * 
 * @param domset The dominating set to check
 * @param N Number of vertices in the graph
 * @param degrees list containing the degree of each vertex
 * @param neighbours the neighbours of each vertex
 * @return true if the given set is secure dominating
 * @return false if the given set is not secure dominating
 */
bool secure_dominates(int* domset, int N, int* degrees, int** neighbours){

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
                    dommed = false;
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


bool connected_dominates(int* domset, int* &dommed, int &domsum, int added, int N, int* degrees, int** neighbours){
    
    return dominates(domset, dommed, domsum, added, N, degrees, neighbours) && connected(domset, N, degrees, neighbours);
}


/**
 * @brief Determine if the given subset is connected.
 * 
 * @param subset The subset to check
 * @param N Number of vertices in the graph
 * @param degrees list containing the degree of each vertex
 * @param neighbours the neighbours of each vertex
 * @return true if the given set is connected
 * @return false if the given set is not connected
 */
bool connected(int* subset, int N, int* degrees, int** neighbours){
    
    int reached[N] = {0};

    int start = -1;

    for (int i = 0; i < N; i++)
    {
        if(subset[i]){
            start = i;
            reached[i] = 1;
            break;
        }
    }
    
    if(start == -1) return false;

    int queue[N] = {0};
    int qsize = 1;
    int qind = 0;
    queue[0] = start;

    while(qind < qsize){
        int curr = queue[qind];
        for (int i = 0; i < degrees[curr]; i++)
        {
            int nb = neighbours[curr][i];
            if(subset[nb] && !reached[nb]){
                reached[nb] = 1;
                queue[qsize++] = nb;
                
            }
        }
        

        qind++;
    }

    for (int i = 0; i < N; i++)
    {
        if(subset[i] && !reached[i]){
            return false;
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
int weight_rand(int N, double* P, double sumP){

    
    double val = (sumP-1e-6)*(rand()/double(RAND_MAX));

    for (int i = 0; i < N; i++)
    {
        if(val < P[i]){
            
            return i;
        }
        val -= P[i];
        
    }

    return 0;
}


int weight_rand_acc(int N, double* P, double sumP, int*& links){
    double val = (sumP)*(rand()/double(RAND_MAX));
    double sum = 1e-6;

    int ind  = links[2*N];
    for(;;){
        // fprintf(stdout, "ind is: %d\n", ind);
        sum += P[ind];

        // fprintf(stdout, "val is: %.10f\n", val);

        if(val <= sum){
            //  fprintf(stdout, "selected: %d\n", ind);
            int prev = links[ind+N];
            int next = links[ind];

            links[prev] = next;
            links[next+N] = prev;

            if(ind == links[2*N]){
                links[2*N] = next;
            }
            return ind;
        }

        ind = links[ind];
        

    }
    
    
    fprintf(stdout, "val is: %.10f, sum is %.10f\n", val, sum);
    fprintf(stdout, "bad times\n");

    return -1;

}
