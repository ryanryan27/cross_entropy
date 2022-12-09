#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cross_entropy.h"


int main(int argc, char* argv[]){

    if(argc == 1){
        print_instructions();
        return 0;
    }
    Params params;
    handle_params(params, argc, argv);

    run_cross_entropy(params);

    return 0;
}

/**
 * @brief Prints instructions on how to run the program from command line
 * 
 */
void print_instructions(){
    fprintf(stdout,
    "   _____                                _                         \n"
    "  / ____|                              | |                        \n"
    " | |     _ __ ___  ___ ___    ___ _ __ | |_ _ __ ___  _ __  _   _ \n"
    " | |    | '__/ _ \\/ __/ __|  / _ \\ '_ \\| __| '__/ _ \\| '_ \\| | | |\n"
    " | |____| | | (_) \\__ \\__ \\ |  __/ | | | |_| | | (_) | |_) | |_| |\n"
    "  \\_____|_|  \\___/|___/___/  \\___|_| |_|\\__|_|  \\___/| .__/ \\__, |\n"
    "                                                     | |     __/ |\n"
    "                                                     |_|    |___/ \n\n\n");

    fprintf(stdout,
    "Runs the cross entropy method on a specified graph to generate one of a few\n"
    "different variants of dominating set (domination, total domination, two-\n"
    "domination, connected domination, or secure domination\n\n");
    fprintf(stdout,
    "This program allows the following options (run the function with no additional\n"
    "arguments to get a help menu):\n\n"
    "  -f <filename>               ---> specify the name of a file containing an\n"
    "     (required)                    edge list with one edge per line\n\n"
    "  -d <char>                   ---> type of dominating set to be generated:\n"
    "     (optional, default d)           d: domination\n"
    "                                     t: total domination\n"
    "                                     2: two-domination\n"
    "                                     c: connected domination\n"
    "                                     s: secure domination\n\n"
    "  -n <int>                    ---> the number of dominating sets  to generate\n"
    "     (optional, default 50)        for each step of the main loop\n\n"
    "  -m <int>                    ---> the number of dominating sets to use to\n"
    "     (optional, default 10)        update the probability distribution\n\n"
    "  -r <int>                    ---> the number of iterations to run without\n"
    "     (optional, default 5)         improvement before the program finishes\n\n"
    "  -R <double>                 ---> number close to zero used to avoid divide\n"
    "     (optional, default 0.1)       by zero errors in score calculations\n\n"
    "  -a <double>                 ---> value between 0 and 1 that specifies how\n"
    "     (optional, default 0.5)       much probability distributions change\n"
    "                                   between iterations\n\n"
    "  -s <int>                    ---> seed for the random function that is\n"
    "     (optional, default 0)         used to select vertices based on the\n"
    "                                   probability distribution\n\n"
    "  -t <double>                 ---> upper limit in runtime in seconds\n"
    "     (optional, default -1)        set to -1 for no limit\n\n"
    "  -o <int>                    ---> output format:\n"
    "     (optional, default 2)          -1: csv format\n"
    "                                     0: no output\n"
    "                                     1: shows domset as a 0-1 vector\n"
    "                                     2: shows domset as list of indices\n"
    "                                     3: same as 1 with final P values below\n\n"
    " -i <int>                     ---> number of seeds to attempt. reports\n"
    "    (optional, default 1)          on the best domset found\n\n"
    );
}


/**
 * @brief Runs the main cross entropy loop with the given param set
 * 
 * @param params set of parameters to be used for cross entropy
 */
void run_cross_entropy(Params params){
    

    //create and initialise the graph to be used for cross entropy, based on the file name specified in params
    Graph graph;
    read_edges(graph, params);
    make_graph(graph, params);

    //create and initialise the data structures used for the cross entropy method
    CEUpdater ce;
    init_updater(ce, graph, params);
    ce.start = clock();

    //run through the cross entropy method for the specified number of seeds
    for (int i = 0; i < params.iterations; i++)
    {
        //if no domset is possible, don't bother trying again
        if (!ce.domset_possible) break;
        
        //set the random seed for this iteration
        srand(params.seed + i);

        //re-initialise the cross entropy data structures
        for (int j = 0;j < graph.N;j++)
        {
            ce.best_domset_this_iteration[j] = 1;
            ce.P[j] = 1.0/double(graph.N);
            ce.Pstar[j] = 0;
        }

        ce.loops_without_change = 0;

        //start the cross entropy main loop
        cross_entropy_main_loop(ce, graph, params);

        //calculate the score of the best domset for this seed
        int sum = 0;
        for (int j = 0; j < graph.N; j++)
        {
            sum += ce.best_domset_this_iteration[j];
        }

        ce.results[i] = sum;

        //keep track of the overall best dominating set for all seeds
        if(ce.results[i] <ce.best){
            ce.best = ce.results[i];
            ce.best_seed = params.seed + i;
            for (int j = 0; j < graph.N; j++)
            {
                ce.best_domset_overall[j] = ce.best_domset_this_iteration[j];
            }
        }

        //stop checking different seeds if we have timed out
        if(ce.timed_out){
            break;
        }
    }

    //check how long the cross entropy method took, for all different seeds
    ce.total_time  = (double)(clock() - ce.start)/CLOCKS_PER_SEC;

    print_output(ce, params, graph);

    destruct_memory(graph, ce, params);

}

/**
 * @brief Run the cross entropy main method using the given graph, parameter set, and initialised updater
 * 
 * @param ce updater that will be updated throughout the loop
 * @param graph graph to generate dominating sets for
 * @param params set of parameters to use for cross entropy
 */
void cross_entropy_main_loop(CEUpdater& ce, Graph graph, Params params){

    while(true){
            
        //generate the domsets based on the current values in P
        //if none are possible, break out of each loop rather then exiting, so that
        //memory can still be freed safely
        for (int j = 0; j < params.n; j++)
        {
            ce.domset_possible = make_domset(ce.domsets[j], graph, ce, params);
            
            if (!ce.domset_possible) break;

            ce.L[j] = calculate_score(graph.N, ce.domsets[j]);
        }

        if (!ce.domset_possible) break;

        //sort the dominating sets based on their calculated scores
        sort_domsets(ce, graph, params);
        
        //if some domset has the best score so far, keep track of it
        if(calculate_score(graph.N, ce.best_domset_this_iteration) > ce.L[0]){
            memcpy(ce.best_domset_this_iteration, ce.domsets[0], graph.N*sizeof(int));
            ce.loops_without_change = 0;
        }

        //if we have had sufficiently many loops without finding a better domset, stop
        if(ce.loops_without_change > params.r){
            break;
        }

        //if we have exceeded the specified time limit, exit
        if(params.timeout > 0 && ((double)(clock()-ce.start)/CLOCKS_PER_SEC) > params.timeout){
            ce.timed_out = true;
            break;
        }

        //calculate the delta value for cross_entropy
        ce.delta = -1*ce.L[0]/log(params.rho);

        //calculate the Pstar value for the set of created domsets
        double psum = 0;
        for (int j = 0; j < graph.N; j++)
        {
            ce.Pstar[j] = calculate_Pstar(j, ce, params);
            psum += ce.Pstar[j];
        }

        //update the values in P based on Pstar and the specified alpha value
        for (int j = 0; j < graph.N; j++)
        {
            ce.P[j] = (1-params.alpha)*ce.P[j] + params.alpha*(ce.Pstar[j]/psum);
        }

        ce.loops_without_change++;
    
    }

}


/**
 * @brief Loads cross entropy parameters based on command line input
 * 
 * @param params paramter struct to be modified
 * @param argc number of specified arguments
 * @param argv character array containing cl arguments
 */
void handle_params(Params& params, int argc, char* argv[]){
    bool has_file = false;
    bool has_out_file = false;
    for (int i = 0; i < argc; i++)
    {
        char* s = argv[i];
        
        if(!strcmp(s, "-f")){
            params.filename = argv[++i];
            has_file = true;
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
            params.timeout = atof(argv[++i]);
        } else if(!strcmp(s, "-o")){
            params.output_types = atof(argv[++i]);

            if(params.output_types == -2 && i+1 < argc && argv[i+1][0] != '-'){

                params.out_file = argv[++i];
                has_out_file = true;
            }

        } else if(!strcmp(s, "-i")){
            params.iterations = atoi(argv[++i]);
        } else if(!strcmp(s, "-d")){

            char* dt = argv[++i];

            params.dom_type_str = dt;

            if(!strcmp(dt, "d")){
                params.dom_type = domination;
            }
            else if(!strcmp(dt, "t")){
                params.dom_type = total;
            }
            else if(!strcmp(dt, "2")){
                params.dom_type = two;
            }
            else if(!strcmp(dt, "s")){
                params.dom_type = secure;
            }
            else if(!strcmp(dt, "c")){
                params.dom_type = connected;
            }
        }
    }
    
    if(!has_file){
        fprintf(stdout, "Please provide a edge list file with -f \"filename\"");
        exit(1);
    }

    if(params.output_types == -2 && !has_out_file){
        fprintf(stdout, "Please specify which file to write to with -o -2 \"filename\"");
        exit(1);
    }
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
    int min_value = INT_MAX;

    while(true){
        int a = 0;
        int b = 0;
        int check = fscanf(file, " %d %d ", &a, &b);
        if(check == EOF) break;

        
        edge_buffer[count] = (int*) malloc(2*sizeof(int));
        edge_buffer[count][0] = a;
        edge_buffer[count++][1] = b;

        if(a < min_value) min_value = a;
        if(b < min_value) min_value = b;
        
    }

    graph.edges = (int**) malloc(count*sizeof(int*));
    for (int i = 0; i < count; i++)
    {
        graph.edges[i] = edge_buffer[i];
        graph.edges[i][0] -= min_value;
        graph.edges[i][1] -= min_value;
    }

    params.label_offset = min_value;
    
    graph.M = count;

    free(edge_buffer);

    
    return count;

}

/**
 * @brief Creates a list of neighbours given a set of edges
 * 
 * @param Graph a graph struct with a specified set of edges.
 */
void make_graph(Graph& graph, Params params){
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

    if(params.dom_type == secure){
        generate_three_apart(graph);
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
    updater.domset_possible = 1;
    updater.timed_out = false;
    updater.best_seed = params.seed;
    updater.best_domset_this_iteration = (int*) malloc(graph.N*sizeof(int));
    memset(updater.best_domset_this_iteration, 1, graph.N*sizeof(int));
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
    updater.best_domset_overall = (int*) malloc(graph.N*sizeof(int));
    memset(updater.best_domset_overall, 1, graph.N*sizeof(int));
}

/**
 * @brief Initialises arrays for the domupdater for the specified paramters and graph
 * 
 * @param du domupdater to initialised
 * @param graph contains array size information
 * @param params contains dom type
 * @param updater contains initial P values
 */
void init_dom_updater(DomUpdater& du, Graph graph, Params params, CEUpdater updater){
    int N = graph.N;

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

    if(params.dom_type == secure){
        du.secure_dommed = (int*)malloc(N*sizeof(int));
        memset(du.secure_dommed, 0, N*sizeof(int));

        du.secure_dom_neighbours = (int**)malloc(N*sizeof(int*));

        for (int i = 0; i < N; i++)
        {
            int size = graph.degrees[i];
            du.secure_dom_neighbours[i] = (int*)malloc(size*sizeof(int));
            memset(du.secure_dom_neighbours[i], 0, size*sizeof(int));
        }

    } else if (params.dom_type == connected){
        du.components = (int*)malloc(N*sizeof(int));
        memset(du.components, -1, N*sizeof(int));
        du.num_components = 0;
        du.subtracting = false;
    }
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

    init_dom_updater(du, graph, params, updater);

    //set up linked list for slightly faster weightrand calculations
    int* links = (int*) malloc((2*N+1)*sizeof(int));
    links[N-1] = -1;
    links[N] = -1;
    links[2*N] = 0;
    
    for(int i = 0; i < N-1; i++){
        links[i] = i+1;
        links[N + i + 1] = i;
    }

    //pre calculate some vertices to add
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
            destruct_memory(du, graph, params);
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


    // set up probabilities for removing vertices
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

    if(params.dom_type == connected){
        du.subtracting = true;
    }

    for (int i = 0; i < N; i++)
    {
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


    destruct_memory(du, graph, params);
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
        double val = pow(exp(1), -updater.L[j]/updater.delta);
        sum_den += val;
        if(updater.domsets[j][i]){
            sum_num += val;
        }
    }
    
    return sum_num/sum_den;
}



/**
 * @brief Print the results of the cross entropy algorithm, based on the specified parameters
 * 
 * @param ce 
 * @param params 
 */
void print_output(CEUpdater ce, Params params, Graph graph){

    FILE* output = stdout;

    if(params.output_types == -2){
        output = fopen(params.out_file, "a");
    }

    if(ce.domset_possible){
        if(params.output_types > 0){
            fprintf(stdout, "Best is %d guards found at seed %d\n", ce.best, ce.best_seed);
            fprintf(stdout, "Time taken: %0.3f \n", ce.total_time);
            fprintf(stdout, "Dominating set: \n");
        }

        for (int i  = 0; i < graph.N; i++)
        {
            if(params.output_types == 1){
                fprintf(stdout, "%d ", ce.best_domset_overall[i]);
            } else if(params.output_types == 2 && ce.best_domset_overall[i]){
                fprintf(stdout, "%d ", i+params.label_offset);
            } else if(params.output_types == 3){
                fprintf(stdout, "%d    ", ce.best_domset_overall[i]);
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

        if(params.output_types < 0){
            fprintf(output, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", params.filename, params.dom_type_str, params.n, params.m, params.r, params.rho, params.alpha, ce.best, ce.total_time);
        }
    } else {
        if(params.output_types < 0){
            fprintf(output, "%s, %s, %d, %d, %d, %f, %.1f, %d, %0.3f\n", params.filename, params.dom_type_str, params.n, params.m, params.r, params.rho, params.alpha, -1, 0.00);
        } else {
            fprintf(stdout, "Unable to dominate graph - probably disconnected\n");
        }
    }

    if(params.output_types == -2){
        fclose(output);
    }
}

/**
 * @brief Free memory allocated for the graph and the cross entropy data
 * 
 * @param graph arrays will be freed
 * @param ce arrays will be freed
 * @param params contains number of domsets to be freed
 */
void destruct_memory(Graph graph, CEUpdater ce, Params params){

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

    
    free(ce.best_domset_this_iteration);
    free(ce.P);
    free(ce.Pstar);
    free(ce.L);
    free(graph.degrees);
    free(graph.edges);
    free(graph.neighbours);
    free(ce.domsets);
    free(ce.results);
    free(ce.best_domset_overall);

    if(params.dom_type == secure){
        for (int i = 0; i < graph.N; i++)
        {
            free(graph.three_aparts[i]);
        }
        free(graph.three_aparts);
        free(graph.degrees3);
        
    }

}

/**
 * @brief Free memory allocated for the domupdater
 * 
 * @param domupdater whose memory is to be freed
 * @param graph contains lengths of arrays to be freed
 * @param params contains dom type
 */
void destruct_memory(DomUpdater du, Graph graph, Params params){

    free(du.dommed);
    free(du.Ptemp);

    if(params.dom_type == secure){
        for (int i = 0; i < graph.N; i++)
        {
            free(du.secure_dom_neighbours[i]);
        }
        free(du.secure_dom_neighbours);
        free(du.secure_dommed);
    } else if(params.dom_type == connected){
        free(du.components);
    }
}


/**
 * @brief Set the domfunc of the domupdater
 * 
 * @param updater updater whose domfunc to set
 * @param params contains the dominating type paramater
 */
void set_domfunc(DomUpdater& updater, Params params){
    if(params.dom_type == total){
        updater.dom_func = &total_dominates;
    }
    else if(params.dom_type == two){
        updater.dom_func = &two_dominates;
    }
    else if(params.dom_type == connected){
        updater.dom_func = &connected_dominates;
    }
    else if(params.dom_type == secure){
        updater.dom_func = &secure_dominates;
    }
    else {
        updater.dom_func = &dominates;
    }
    
}


/**
 * @brief Determine if the given dominating set dominates the graph defined by the list of neighbours.
 * 
 * @param du the domupdated used to keep track of domination info
 * @param added the vertex where a guard was changed
 * @param domset The dominating set to check
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
            } 
        }
    }

    return du.domsum == graph.N;

}

/**
 * @brief Determine if the given dominating set total dominates the graph defined by the list of neighbours.
 * 
 * @param du the domupdated used to keep track of domination info
 * @param added the vertex where a guard was changed
 * @param domset The dominating set to check
 * @param graph the graph to check domination over
 */
bool total_dominates(DomUpdater& du, int added, int* domset, Graph graph){

    if(domset[added]){

        for (int j = 0; j < graph.degrees[added]; j++)
        {
            if(!du.dommed[graph.neighbours[added][j]]) 
            {
                du.domsum++;
                du.dommed[graph.neighbours[added][j]] = 1;
            }
        }

    } else {

        for (int i = 0; i < graph.degrees[added]; i++)
        {
            int neighbour = graph.neighbours[added][i];

            bool undommed = true;
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
            } 
        }
    }

    return du.domsum == graph.N;

}


/**
 * @brief Determine if the given dominating set two-dominates the graph defined by the list of neighbours.
 * 
 * @param du the domupdated used to keep track of domination info
 * @param added the vertex where a guard was changed
 * @param domset The dominating set to check
 * @param graph the graph to check domination over
 */
bool two_dominates(DomUpdater& du, int added, int* domset, Graph graph){

    if(domset[added]){
        if(du.dommed[added] < 2) du.domsum++;
        du.dommed[added] = 2;
        for (int j = 0; j < graph.degrees[added]; j++)
        {
            if(!domset[graph.neighbours[added][j]]) 
            {
                du.dommed[graph.neighbours[added][j]]++;

                if(du.dommed[graph.neighbours[added][j]] == 2) du.domsum++;
                
            }
        }

    } else {

        int domset_nbrs = 0;
        for (int j = 0; j < graph.degrees[added]; j++)
        {
            int nbr = graph.neighbours[added][j];
            if(domset[nbr]){
                domset_nbrs++;
            }
        }
        du.dommed[added] = domset_nbrs;
        if(domset_nbrs < 2){
            du.domsum--;
            return false;
        } 


        for (int i = 0; i < graph.degrees[added]; i++)
        {
            int neighbour = graph.neighbours[added][i];

            if(domset[neighbour]) continue;

            domset_nbrs = 0;
            for (int j = 0; j < graph.degrees[neighbour]; j++)
            {
                int nbr = graph.neighbours[neighbour][j];
                if(domset[nbr]){
                    domset_nbrs++;
                }
            }
            
            du.dommed[neighbour] = domset_nbrs;
             
            if(domset_nbrs < 2){
                du.domsum--;
                return false;
            } 
        }
    }

    return du.domsum == graph.N;

}

/**
 * @brief Determine if the given dominating set two-dominates the graph defined by the list of neighbours.
 * 
 * @param du the domupdated used to keep track of domination info
 * @param added the vertex where a guard was changed
 * @param domset The dominating set to check
 * @param graph the graph to check domination over
 */
bool secure_dominates(DomUpdater& du, int added, int* domset, Graph graph){

    //if the vertex in question has been added to the domset
    if(domset[added]){

        //update which vertices are dominated by added
        du.dommed[added]++;

        for (int i = 0; i < graph.degrees[added]; i++)
        {
            du.dommed[graph.neighbours[added][i]]++;
        }
        
        //check which vertices within three edges of the new guard are now dominated
        for (int i = 0; i < graph.degrees3[added]; i++)
        {
            int vertex_within_three = graph.three_aparts[added][i];

            //only consider guards within 3
            if(!domset[vertex_within_three]) continue;

            //check the neighbours of these guards
            for (int j = 0; j < graph.degrees[vertex_within_three]; j++)
            {
                int within_three_neighbour = graph.neighbours[vertex_within_three][j];

                //if was secure dominated before, it still will be
                if(du.secure_dom_neighbours[vertex_within_three][j]) continue;

                //if it is now dominated, update necessary info
                if(can_secure_dom(vertex_within_three, within_three_neighbour, domset, graph, du)){
                    du.secure_dom_neighbours[vertex_within_three][j] = 1;
                    du.secure_dommed[within_three_neighbour]++;
                }
            }
        }
    }
    //if the vertex in question has been removed from the domset...
    else{

        //update which vertices are dominated by added
        du.dommed[added]--;

        for (int i = 0; i < graph.degrees[added]; i++)
        {
            du.dommed[graph.neighbours[added][i]]--;
        }

        //check which of the vertices within 3 edges are still secure dominated
        for (int i = 0; i < graph.degrees3[added]; i++)
        {
            int vertex_within_three = graph.three_aparts[added][i];

            //only interested in guards and the vertex that was removed
            if(vertex_within_three != added && !domset[vertex_within_three]) continue;

            //check the neighbours of the guards or removed guard
            for (int j = 0; j < graph.degrees[vertex_within_three]; j++)
            {
                int within_three_neighbour = graph.neighbours[vertex_within_three][j];

                //if wasnt secure dommed before, still wont be
                if(!du.secure_dom_neighbours[vertex_within_three][j]) continue;

                //if it is no longer secure dommed, update necessary info
                if(!can_secure_dom(vertex_within_three, within_three_neighbour, domset, graph, du)){
                    du.secure_dom_neighbours[vertex_within_three][j] = 0;
                    du.secure_dommed[within_three_neighbour]--;
                }
            }
        }

    }

    //check that everything is either secure dominated by a neighbour or by itself
    for (int i = 0; i < graph.N; i++)
    {
        if(domset[i] == 0 && du.secure_dommed[i] == 0) return false;
    }

    return true;
    
}


/**
 * @brief Determine if the given dominating set connected-dominates the graph defined by the list of neighbours.
 * 
 * @param du the domupdated used to keep track of domination info
 * @param added the vertex where a guard was changed
 * @param domset The dominating set to check
 * @param graph the graph to check domination over
 */
bool connected_dominates(DomUpdater& du, int added, int* domset, Graph graph){
    bool doms = dominates(du, added, domset, graph);
    bool conn;

    if(!du.subtracting){
        conn = is_connected(du, added, domset, graph);
    } else if(domset[added]){
        conn = 1;
    } else {
        conn = is_connected(domset, graph);
    }

    // if(!du.subtracting && conn != is_connected(domset, graph)){
    //     fprintf(stdout, "not the same boss\n");

    //     for (int i = 0; i < graph.N; i++)
    //     {
    //         fprintf(stdout, "%d ", du.components[i]+1);
    //     }
    //     fprintf(stdout, "\n");
        

    //     exit(0);
    // }

    return doms && conn;
}

bool is_connected(DomUpdater& du, int added, int* domset, Graph graph){

    int degree = graph.degrees[added];

    int num_adjacent = 0;
    int* others = (int*)malloc(degree*sizeof(int));

    for (int i = 0; i < degree; i++)
    {
        int v2 = graph.neighbours[added][i];
        int component_v2 = du.components[v2];

        if(!domset[v2]) continue;

        bool in_already = false;
        for (int j = 0; j < num_adjacent; j++)
        {
            if(others[j] == component_v2){
                in_already = true;
                break;
            }
        }

        if(in_already) continue;

        others[num_adjacent] = component_v2;
        num_adjacent++;
        
    }
    
    if(num_adjacent == 0){
        du.components[added] = added;
        du.num_components++;
    } else if(num_adjacent == 1){
        du.components[added] = others[0];
    } else {
        for (int i = 0; i < graph.N; i++)
        {
            int component_i = du.components[i];
            for (int j = 1; j < num_adjacent; j++)
            {
                if(component_i == others[j]){
                    du.components[i] = others[0];
                    break;
                }
            }    
        }

        du.components[added] = others[0];
        du.num_components += 1 - num_adjacent;

    }


    free(others);

    //fprintf(stdout, "num_components: %d, num_adjacent: %d\n", du.num_components, num_adjacent);

    return du.num_components == 1;
}

/**
 * @brief Determine if the given subset is connected in the given graph;
 * 
 * @param subset The subset of vertices whose induced subgraph will be checked
 * @param graph The graph to check connection under
 */
bool is_connected(int* subset, Graph graph){
    
    //just do a breadth first search only considering vertices from the subset
    int reached[graph.N] = {0};

    int start = -1;

    for (int i = 0; i < graph.N; i++)
    {
        if(subset[i]){
            start = i;
            reached[i] = 1;
            break;
        }
    }
    
    if(start == -1) return false;

    int queue[graph.N] = {0};
    int qsize = 1;
    int qind = 0;
    queue[0] = start;
    int num_reached = 0;
    while(qind < qsize){
        int curr = queue[qind];
        for (int i = 0; i < graph.degrees[curr]; i++)
        {
            int nb = graph.neighbours[curr][i];
            if(subset[nb] && !reached[nb]){
                reached[nb] = 1;
                queue[qsize++] = nb;
                num_reached++;
                
            }
        }
        

        qind++;
    }

    //fprintf(stdout, "num_reached from start at %d: %d\n", start, num_reached);
    
    for (int i = 0; i < graph.N; i++)
    {
        if(subset[i] && !reached[i]){
            return false;
        }
    }

    return true;
}

/**
 * @brief Select an index between 0 and N based on the relative weight of values in P with that index.
 * 
 * @param N max number of the index
 * @param P list of weights for each index, higher means more likely to be chosen
 * @param sumP sum of all elements in P
 * @param links linked list that keeps track of which elements of P have already been selected
 * @return the chosen index between 0 and N
 */
int weight_rand_acc(int N, double* P, double sumP, int*& links){

    //generate a random value between 0 and the sum of P
    double val = (sumP)*(rand()/double(RAND_MAX));
    
    //offset for very small values
    double sum = 1e-6;

    //iterate through the linked list until the cumulative values of P exceed the randomly chosen value
    int ind  = links[2*N];
    for(;;){
        sum += P[ind];

        //once we've chosen one, update the linked list to skip that index next time
        if(val <= sum){
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

}


/**
 * @brief Checks if a guard at a vertex can securely dominate the target vertex
 * for a given domset, graph and domupdater
 * 
 * @param guard vertex to check
 * @param target target vertex to check
 * @param domset domset used to check
 * @param graph graph to check for
 * @param du updater containing information about the current status of other vertices relevant to secure domiantion
 * @return true if target can be securely dominated by guard
 * @return false otherwise
 */
bool can_secure_dom(int guard, int target, int* domset, Graph graph, DomUpdater du){
    if(!domset[guard]) return false;

    for (int i = 0; i < graph.degrees[guard]; i++)
    {
        int guard_neighbour = graph.neighbours[guard][i];

        if(guard_neighbour == target) continue;

        if(domset[guard_neighbour]) continue;

        if(du.dommed[guard_neighbour] > 1) continue;

        bool is_neighbour = false;
        for (int j = 0; j < graph.degrees[target]; j++)
        {
            if(guard_neighbour == graph.neighbours[target][j]) is_neighbour = true;
        }
        
        if(is_neighbour) continue;

        return false;
    }

    return true;
    
}


/**
 * @brief creates an array of arrays that specify which vertices are within 3 edges of a given vertex
 *  this includes the original vertex
 * 
 * @param graph graph to create the array for
 */
void generate_three_apart(Graph& graph){
    int N = graph.N;

    graph.degrees3 = (int*)malloc(N*sizeof(int));
    graph.three_aparts = (int**)malloc(N*sizeof(int*));

    int one_away[N];
    int two_away[N];
    int three_away[N];

    for (int i = 0; i < N; i++)
    {
        memset(one_away, 0, N*sizeof(int));
        memset(two_away, 0, N*sizeof(int));
        memset(three_away, 0, N*sizeof(int));

        int within_three_sum = 1;

        //get which vertices are one away
        for (int j = 0; j < graph.degrees[i]; j++){
            one_away[graph.neighbours[i][j]] = 1;
            within_three_sum++;
        }

        //get which vertice are exactly two away
        for (int j = 0; j < graph.degrees[i]; j++)
        {
            int neighbour = graph.neighbours[i][j];
            for (int k = 0; k < graph.degrees[neighbour]; k++)
            {
                int second_neighbour = graph.neighbours[neighbour][k];

                //accounts for triangles and not adding self again
                if(one_away[second_neighbour] || two_away[second_neighbour] || second_neighbour == i) continue;

                two_away[second_neighbour] = 1;
                within_three_sum++;

            }
            
        }

        //get vertices which are exactly three away
        for (int j = 0; j < N; j++)
        {
            if(two_away[j] == 0) continue;

            for (int k = 0; k < graph.degrees[j]; k++)
            {
                int third_neighbour = graph.neighbours[j][k];

                //accounts for C4s
                if(one_away[third_neighbour] || two_away[third_neighbour] || three_away[third_neighbour] || third_neighbour == i) continue;

                three_away[third_neighbour] = 1;
                within_three_sum++;

            }
            
        }

        graph.degrees3[i] = within_three_sum;

        graph.three_aparts[i] = (int*) malloc(within_three_sum*sizeof(int));

        int counter = 0;

        for (int j = 0; j < N; j++)
        {
            if(j == i || one_away[j] || two_away[j] || three_away[j]){
                graph.three_aparts[i][counter] = j;
                counter++;
            }
        }
        
        if(counter != within_three_sum){
            fprintf(stdout, "three apart not match: counter=%d, degree=%d\n", counter, within_three_sum);
            exit(1);
        }
        
        
    }
    

}