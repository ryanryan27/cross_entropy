#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <math.h>
#include "cross_entropy.h"

static const double e = 2.71828;


int main(int argc, char* argv[]){

    //Graph Variables
    std::string filename;
    int N;
    int max_degree;

    int* degrees;
    int** neighbours;
    int** edges; //buffer for reading in edges.

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
    
    bool has_file = 0;

    for (int i = 0; i < argc; i++)
    {
        std::string s(argv[i]);
        if(s=="-f"){
            filename = argv[++i];
            has_file = 1;

            if(i + 1 < argc && argv[i+1][0] != '-'){
                label_offset = atoi(argv[++i]);
            }
        } else if(s == "-n"){
            n = atoi(argv[++i]);
        } else if(s == "-m"){
            m = atoi(argv[++i]);
        } else if(s == "-r"){
            r = atoi(argv[++i]);
        } else if(s == "-R"){
            rho = atof(argv[++i]);
        } else if(s == "-a"){
            alpha = atof(argv[++i]);
        }
    }
    
    if(!has_file){
        std::cout << "Please provide a edge list file with -f \"filename\"" << std::endl;
        return 1;
    }

    int edge_count = read_edges(edges, filename, label_offset);

    make_graph(N, degrees, neighbours, edge_count, edges);
    


    dombest = new int[N];
    P = new double[N];
    Pstar = new double[N];
    for (int i = 0;i < N;i++)
    {
        dombest[i] = 1;
        P[i] = 1.0/double(N);
        Pstar[i] = 0;
    }

    int t = 0;
    while(true){

        domsets = new int*[n];
        L =  new double[n];
        for (int i = 0; i < n; i++)
        {
            make_domset(domsets[i], N, degrees, neighbours, P);
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
            
    std::cout << "Best is " << sum << " guards" << std::endl;


    return 0;
}


int read_edges(int**& edges, std::string filename, int label_offset){
    std::ifstream file(filename.c_str());

    if(!file){
        std::cout << "File does not exist" << std::endl;
        return -1;
    }

    std::string line;

    int** edge_buffer = new int*[10000];

    int count = 0;
    while (std::getline(file, line)){
        int space = line.find(" ");
        edge_buffer[count] = new int[2];
        edge_buffer[count][0] = atoi(line.substr(0,space).c_str()) - label_offset;
        edge_buffer[count++][1] = atoi(line.substr(space+1,line.length()-1).c_str()) - label_offset;
    }

    edges = new int*[count];
    for (int i = 0; i < count; i++)
    {
        edges[i] = new int[2];
        edges[i][0] = edge_buffer[i][0];
        edges[i][1] = edge_buffer[i][1];
    }
    

    
    return count;

}

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

    degrees = new int[N];
    memset(degrees, 0, N*sizeof(int));
    

    for (int i = 0; i < edge_count; i++)
    {
        int a = edges[i][0];
        int b = edges[i][1];

        degrees[a]++;
        degrees[b]++;
    }


    

    neighbours = new int*[N];
    for (int i = 0; i < N; i++)
    {
        neighbours[i] = new int[degrees[i]];
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

void make_domset(int* &domset, int N, int* degrees, int** neighbours, double* P){
    domset = new int[N];
    memset(domset, 0, N*sizeof(int));

    double Ptemp[N];
    std::copy(P, P+N, Ptemp);

    while(!dominates(domset, N, degrees, neighbours)){
        int ind = weight_rand(N, Ptemp);
        domset[ind] = 1;
        Ptemp[ind] = 0;
    }

    for (int i = 0; i < N; i++){
        if(domset[i]){
            domset[i] = 0;
            if(!dominates(domset, N, degrees, neighbours)){
                domset[i] = 1;
            }
        }
    }


}

double calculate_score(int N, int* domset){
    int sum = 0;

    for (int i = 0; i < N; i++)
    {
        sum += domset[i];
    }
    
    return sum;
}

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

int compare_scores( const void* a, const void* b)
{
    double L_a = ((double*)a)[1];
    double L_b = ((double*)b)[1];

     if ( L_a == L_b ) return 0;
     else if ( L_a < L_b ) return -1;
     else return 1;
}

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