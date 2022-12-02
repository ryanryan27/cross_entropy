#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cross_entropy.h"

static const double e = 2.71828;


int main(int argc, char* argv[]){


    Graph graph;

    Params params;

    CEUpdater ce;

    double timeout = -1;
    
    bool has_file = 0;

    for (int i = 0; i < argc; i++)
    {
        char* s = argv[i];
        
        if(!strcmp(s, "-f")){
            params.filename = argv[++i];
            has_file = 1;

            if(i + 1 < argc && argv[i+1][0] != '-'){
                params.label_offset = atoi(argv[++i]);
            }
        } else if(!strcmp(s, "-n")){
            params.n = atoi(argv[++i]);
        } else if(!strcmp(s, "-m")){
            params.m = atoi(argv[++i]);
        } else if(!strcmp(s, "-r")){
            params.r = atoi(argv[++i]);
        } else if(!strcmp(s, "-R")){
            params.rho = atof(argv[++i]);
        } else if(!strcmp(s, "-a")){
            params.alpha = atof(argv[++i]);
        } else if(!strcmp(s, "-s")){
            params.seed = atof(argv[++i]);
        } else if(!strcmp(s, "-t")){
            timeout = atof(argv[++i]);
        } else if(!strcmp(s, "-o")){
            params.output_types = atof(argv[++i]);
        } else if(!strcmp(s, "-i")){
            params.iterations = atoi(argv[++i]);
        } else if(!strcmp(s, "-d")){
            params.dom_type = argv[++i];
        }
    }
    
    if(!has_file){
        fprintf(stdout, "Please provide a edge list file with -f \"filename\"");
        return 1;
    }

    clock_t start = clock();

    graph.M = read_edges(graph, params);


    make_graph(graph);


    init_updater(ce, graph, params);
    
    bool timed_out = false;
    int its_completed = 0;
    int domset_possible = 1;

    
    for (int i = 0; i < params.iterations; i++)
    {
        if (!domset_possible) break;
        
        srand(params.seed + i);

        for (int j = 0;j < graph.N;j++)
        {
            ce.dombest[j] = 1;
            ce.P[j] = 1.0/double(graph.N);
            ce.Pstar[j] = 0;
        }

        int t = 0;

        while(true){
            
            for (int j = 0; j < params.n; j++)
            {
                domset_possible = make_domset(ce.domsets[j], graph, ce, params);
                
                if (!domset_possible) break;

                ce.L[j] = calculate_score(graph.N, ce.domsets[j]);
            }



            if (!domset_possible) break;

            sort_domsets(ce, graph, params);
            
            if(calculate_score(graph.N, ce.dombest) > ce.L[0]){
                memcpy(ce.dombest, ce.domsets[0], graph.N*sizeof(int));
                t = 0;
            }

            if(t > params.r){
                break;
            }

            if(timeout > 0 && ((double)(clock()-start)/CLOCKS_PER_SEC) > timeout){
                timed_out = true;
                break;
            }


            ce.delta = -1*ce.L[0]/log(params.rho);

            double psum = 0;
            for (int j = 0; j < graph.N; j++)
            {
                ce.Pstar[j] = calculate_Pstar(j, ce, params);
                psum += ce.Pstar[j];
            }

            for (int j = 0; j < graph.N; j++)
            {
                ce.P[j] = (1-params.alpha)*ce.P[j] + params.alpha*(ce.Pstar[j]/psum);
            }

            t++;
        
        }


        int sum = 0;
        for (int j = 0; j < graph.N; j++)
        {
            sum += ce.dombest[j];
        }

        ce.results[i] = sum;

        if(ce.results[i] <ce.best){
            ce.best = ce.results[i];
            for (int j = 0; j < graph.N; j++)
            {
                ce.best_domset[j] = ce.dombest[j];
            }
            
        }

        its_completed = i+1;
        if(timed_out){
            
            break;
        }
    }



    double total_time  = (double)(clock() - start)/CLOCKS_PER_SEC;
    if(domset_possible){
        if(params.output_types > 0){
            fprintf(stdout, "Best is %d guards\n", ce.best);
            fprintf(stdout, "Time taken: %0.3f \n", total_time);
            fprintf(stdout, "Dominating set: \n");
        }

        for (int i  = 0; i < graph.N; i++)
        {
            if(params.output_types == 1){
                fprintf(stdout, "%d ", ce.best_domset[i]);
            } else if(params.output_types == 2 && ce.best_domset[i]){
                fprintf(stdout, "%d ", i+params.label_offset);
            } else if(params.output_types == 3){
                fprintf(stdout, "%d    ", ce.best_domset[i]);
            }
            
        }
        if(params.output_types > 0){
            fprintf(stdout, "\n");
        }
        if(params.output_types == 3){
            for (int i  = 0; i < graph.N; i++)
            {
                fprintf(stdout, "%.2f ", ce.P[i]);
            }
        
            fprintf(stdout, "\n");
        }

        if(params.output_types == -1){
            fprintf(stdout, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", params.filename, params.dom_type, params.n, params.m, params.r, params.rho, params.alpha, ce.best, total_time);
        }
    } else {
        if(params.output_types == -1){
            fprintf(stdout, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", params.filename, params.dom_type, params.n, params.m, params.r, params.rho, params.alpha, -1, 0.00);
        } else {
            fprintf(stdout, "Unable to dominate graph - probably disconnected\n");
        }
    }

    
    for (int i = 0; i < graph.M; i++)
    {
        free(graph.edges[i]);
    }
    for(int i = 0; i < graph.N; i++)
    {
        free(graph.neighbours[i]);
    }
    for (int i = 0; i < params.n; i++)
    {
        free(ce.domsets[i]);
    }

    
    free(ce.dombest);
    free(ce.P);
    free(ce.Pstar);
    free(ce.L);
    free(graph.degrees);
    free(graph.edges);
    free(graph.neighbours);
    free(ce.domsets);
    free(ce.results);
    free(ce.best_domset);


    return 0;
}

/**
 * @brief Reads a list of edges from a given file. 
 * 
 * @param graph a graph struct whose edges will be updated to contain those listed in the specified file
 * @param Params a set of parameters that contains the filename of the edge list file
 */
int read_edges(Graph& graph, Params params){
    FILE* file = fopen(params.filename, "r");

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
        edge_buffer[count][0] = a - params.label_offset;
        edge_buffer[count++][1] = b - params.label_offset;
        
    }

    graph.edges = (int**) malloc(count*sizeof(int*));
    for (int i = 0; i < count; i++)
    {
        graph.edges[i] = edge_buffer[i];
    }
    
    free(edge_buffer);

    
    return count;

}

/**
 * @brief Creates a list of neighbours given a set of edges
 * 
 * @param Graph a graph struct with a specified set of edges.
 */
void make_graph(Graph& graph){
    graph.N = -1;
    for (int i = 0; i < graph.M; i++)
    {
        if(graph.edges[i][0] > graph.N){
            graph.N = graph.edges[i][0];
        }
        if(graph.edges[i][1] > graph.N){
            graph.N = graph.edges[i][1];
        }
    }
    graph.N++;

    graph.degrees = (int*) malloc(graph.N*sizeof(int));
    memset(graph.degrees, 0, graph.N*sizeof(int));
    
    for (int i = 0; i < graph.M; i++)
    {
        int a = graph.edges[i][0];
        int b = graph.edges[i][1];

        graph.degrees[a]++;
        graph.degrees[b]++;
    }

    graph.neighbours = (int**) malloc(graph.N*sizeof(int*));
    for (int i = 0; i < graph.N; i++)
    {
        graph.neighbours[i] = (int*) malloc(graph.degrees[i]*sizeof(int));
    }
    
    int counter[graph.N] = {0};
    for (int i = 0; i < graph.M; i++)
    {

        int a = graph.edges[i][0];
        int b = graph.edges[i][1];


        graph.neighbours[a][counter[a]++] = b;
        graph.neighbours[b][counter[b]++] = a;

    }
    
}

/**
 * @brief Initialise arrays for the ce update tracking
 * 
 * @param updater will have its arrays initialised
 * @param graph data will be used to set array sizes
 * @param params data will be used to set array sizes
 */
void init_updater(CEUpdater& updater, Graph graph, Params params){
    updater.dombest = (int*) malloc(graph.N*sizeof(int));
    memset(updater.dombest, 1, graph.N*sizeof(int));
    updater.P = (double*) malloc(graph.N*sizeof(double));
    updater.Pstar = (double*) malloc(graph.N*sizeof(double));

    updater.results = (int*) malloc(params.iterations*sizeof(int));

    updater.domsets = (int**) malloc(params.n*sizeof(int*));
    updater.L =  (double*) malloc(params.n*sizeof(double));

    for (int i = 0; i < params.n; i++)
    {
        updater.domsets[i] = (int*) malloc(graph.N*sizeof(int));
    }

    updater.best = graph.N;
    updater.best_domset = (int*) malloc(graph.N*sizeof(int));
    memset(updater.best_domset, 1, graph.N*sizeof(int));
}

/**
 * @brief Fills a given array with a dominating set
 * 
 * @param domset array that will contain a domset
 * @param graph the graph to make a domset for
 * @param updater contains probability distribution for vertex selection
 * @param params contains cross entropy parameters
 */
int make_domset(int* &domset, Graph graph, CEUpdater updater, Params params){
    DomUpdater du;

    set_domfunc(du, params);
    
    int N = graph.N;
    memset(domset, 0, N*sizeof(int));

    du.dommed = (int*) malloc(N*sizeof(int));
    memset(du.dommed, 0, N*sizeof(int));

    du.domsum = 0;

    du.Ptemp = (double*) malloc(N*sizeof(double));
    memcpy(du.Ptemp, updater.P, N*sizeof(double));

    du.sumP = 0;
    

    for (int i = 0; i < N; i++)
    {
        du.Ptemp[i] += 1;
        du.sumP += du.Ptemp[i];
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

    int pre_calc = 1;//floor(sqrt(N/2));

    int choices[pre_calc];

    for (int i = 0; i < pre_calc; i++)
    {
        choices[i] = weight_rand_acc(N, du.Ptemp, du.sumP, links);
        du.sumP -= du.Ptemp[choices[i]];
        du.Ptemp[choices[i]] = 0;
    }
    

    int choice_num = 0;


    int domcount = 0;
    int ind;
    do {
        if(domcount == N){
            free(du.dommed);
            free(du.Ptemp);
            free(links);
            return 0;    
        }

        // ind = weight_rand(N, Ptemp, sumP);

        if(choice_num >= pre_calc){
            for (int i = 0; i < pre_calc; i++)
            {
                int indx = weight_rand_acc(N, du.Ptemp, du.sumP, links);

                choices[i] = indx;
                if(indx == -1) break;

                du.sumP -= du.Ptemp[indx];
                du.Ptemp[indx] = 0;
            }
            choice_num = 0;
        }

        ind = choices[choice_num++];

        domset[ind] = 1;
        domcount++;
        

    } while(!du.dom_func(du, ind, domset, graph));


    for (int i = 0; i < N; i++){

        du.Ptemp[i] = 1.0 - updater.P[i];

        if(du.Ptemp[i] < 0 ){
            fprintf(stdout, "oh no \n");
            exit(0);
        }

    }


    du.sumP = 0;

    for (int i = 0; i < N; i++)
    {
        du.sumP += du.Ptemp[i];
    }

    links[N-1] = -1;
    links[N] = -1;
    links[2*N] = 0;
    
    for(int i = 0; i < N-1; i++){
        links[i] = i+1;
        links[N + i + 1] = i;
    }


    for (int i = 0; i < pre_calc; i++)
    {
        choices[i] = weight_rand_acc(N, du.Ptemp, du.sumP, links);
        du.sumP -= du.Ptemp[choices[i]];
        du.Ptemp[choices[i]] = 0;
    }
    

    choice_num = 0;

    //  fprintf(stdout, "removing \n");

    for (int i = 0; i < N; i++)
    {
        // int ind = weight_rand(N, Ptemp, sumP);
        if(choice_num >= pre_calc){
            for (int i = 0; i < pre_calc; i++)
            {
                int indx = weight_rand_acc(N, du.Ptemp, du.sumP, links);

                choices[i] = indx;
                if(indx == -1) break;

                du.sumP -= du.Ptemp[indx];
                du.Ptemp[indx] = 0;
            }
            choice_num = 0;
        }

        int ind = choices[choice_num++];

        if(domset[ind]){
            domset[ind] = 0;
            if(!du.dom_func(du, ind, domset, graph)){
                domset[ind] = 1;
                du.dom_func(du, ind, domset, graph);
            }
        }

    }

    //fprintf(stdout, "here \n");

    free(du.dommed);
    free(du.Ptemp);
    free(links);
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
 * @param updater will have its L and Domsets arrays sorted based on the values in L
 */
void sort_domsets(CEUpdater& updater, Graph graph, Params params){
    double to_sort[params.n][2];
    for (int i = 0; i < params.n; i++)
    {
        to_sort[i][0] = i;
        to_sort[i][1] = updater.L[i];
    }
    
    qsort(to_sort, params.n, sizeof(*to_sort), compare_scores);
    
    double Ltemp[params.n];
    int domset_temp[params.n][graph.N];

    for (int i = 0; i < params.n; i++)
    {
        Ltemp[i] = updater.L[(int)to_sort[i][0]];
        for (int j = 0; j < graph.N; j++)
        {
            domset_temp[i][j] = updater.domsets[(int)to_sort[i][0]][j];
        }
        
    }

    for (int i = 0; i < params.n; i++)
    {
        updater.L[i] = Ltemp[i];
        for (int j = 0; j < graph.N; j++)
        {
            updater.domsets[i][j] = domset_temp[i][j];
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
 * @param updater contains the current cross entropy values
 * @param params contains the cross entropy parameters
 */
double calculate_Pstar(int i, CEUpdater updater, Params params){

    double sum_num = 0;
    double sum_den = 0;
    
    for (int j = 0; j < params.m; j++)
    {
        double val = pow(e, -updater.L[j]/updater.delta);
        sum_den += val;
        if(updater.domsets[j][i]){
            sum_num += val;
        }
    }
    
    return sum_num/sum_den;
}

/**
 * @brief Set the domfunc of the domupdater
 * 
 * @param updater updater whose domfunc to set
 * @param params contains the dominating type paramater
 */
void set_domfunc(DomUpdater& updater, Params params){
    updater.dom_func = &dominates;
}


/**
 * @brief Determine if the given dominating set dominates the graph defined by the list of neighbours.
 * 
 * @param domset The dominating set to check
 * @param added the vertex where a guard was changed
 * @param du the domupdated used to keep track of domination info
 * @param graph the graph to check domination over
 */
bool dominates(DomUpdater& du, int added, int* domset, Graph graph){

    if(domset[added]){
        if(!du.dommed[added]) du.domsum++;
        du.dommed[added] = 1;
        for (int j = 0; j < graph.degrees[added]; j++)
        {
            if(!du.dommed[graph.neighbours[added][j]]) 
            {
                du.domsum++;
                du.dommed[graph.neighbours[added][j]] = 1;
            }
        }

    } else {

        bool undommed = true;
        for (int j = 0; j < graph.degrees[added]; j++)
        {
            int nbr = graph.neighbours[added][j];
            if(domset[nbr]){
                undommed = false;
                break;
            }
        }

        if(undommed){
            du.dommed[added] = 0;
            du.domsum--;
            //return false;
        } 

        for (int i = 0; i < graph.degrees[added]; i++)
        {
            int neighbour = graph.neighbours[added][i];

            if(domset[neighbour]) continue;

            undommed = true;
            for (int j = 0; j < graph.degrees[neighbour]; j++)
            {
                int nbr = graph.neighbours[neighbour][j];
                if(domset[nbr]){
                    undommed = false;
                    break;
                }
            }

            if(undommed){
                du.dommed[neighbour] = 0;
                du.domsum--;
                //return false;
            } 
        }
    }

    return du.domsum == graph.N;

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
    
    return false;//dominates(domset, dommed, domsum, added, N, degrees, neighbours) && connected(domset, N, degrees, neighbours);
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
        //fprintf(stdout, "ind is: %d\n", ind);
        sum += P[ind];
        // fprintf(stdout, "val is: %.10f\n", val);

        if(val <= sum){
            //  fprintf(stdout, "selected: %d\n", ind);
            int prev = links[ind+N];
            int next = links[ind];
            if(prev >= 0){
                links[prev] = next;
            }
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
