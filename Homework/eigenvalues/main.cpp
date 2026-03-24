#include"matrix.hpp" 
#include<vector> 
#include<iostream>
#include<random>
#include<cmath>
#include<fstream>
#include<thread>


pp::matrix random_symmetric_matrix(int n) {
    // Random device for seeding
    std::random_device rd;

    // Mersenne Twister engine
    std::mt19937 gen(rd());

    // Uniform distribution between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    pp::matrix A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A.cols[j][i] = dist(gen);
        }
    }
    // make it symmetric: A = (A + A^T) / 2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            A[i][j] = A[j][i] = 0.5 * (A[i][j] + A[j][i]);
        }
    }
    return A;
}

pp::vector random_vector(int n) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    pp::vector b(n);

    for (int j=0; j<n; j++) {
        b[j] = dist(gen);
    }
    return b;
}

// linspace like in python to make life easier
std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result;

    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}

// jacobi struct
struct jacobi{
    int timesJ(pp::matrix& A, int p, int q, double theta) {  //pass A as reference else while loop gets stuck
       	double c=cos(theta),s=sin(theta);
	    for(int i=0;i<A.size1();i++){
            double aip=A[i][p],aiq=A[i][q];
            A[i][p]=c*aip-s*aiq;
            A[i][q]=s*aip+c*aiq;
		    }
        return 0;
    }
    int Jtimes(pp::matrix& A, int p, int q, double theta) {  // pass A as reference else while loop gets stuck
       	double c=cos(theta),s=sin(theta);
        for(int j=0;j<A.size1();j++){
            double apj=A[p][j],aqj=A[q][j];
            A[p][j]= c*apj+s*aqj;
            A[q][j]=-s*apj+c*aqj;
            }
        return 0;
    }

    std::pair<pp::vector, pp::matrix> cyclic(pp::matrix M) {
        pp::matrix A=M.copy();
        pp::matrix V=M.copy();
        V.setid();
        pp::vector w(M.size1());
        /* run Jacobi rotations on A and update V */
        /* copy diagonal elements into w */
        bool changed;
        do{
            changed=false;
            int n = M.size1();
            for(int p=0;p<n-1;p++)
            for(int q=p+1;q<n;q++){
                double apq=A[p][q], app=A[p][p], aqq=A[q][q];
                double theta=0.5*std::atan2(2*apq,aqq-app);
                double c=std::cos(theta),s=std::sin(theta);
                double new_app=c*c*app-2*s*c*apq+s*s*aqq;
                double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
                if(new_app!=app || new_aqq!=aqq) {// do rotation
                    changed=true;
                    timesJ(A,p,q, theta); // A←A*J 
                    Jtimes(A,p,q, -theta); // A←JT*A 
                    timesJ(V,p,q, theta); // V←V*J
                    }
            }
        }while(changed);

        for (int i = 0; i < A.size1(); i++) {
            w[i] = A[i][i];
        }
        return {w, V};
        }
};

//  result struct for prallel processing
struct timing_result {
    int n;
    double time;
};

// function, which times jacobi struct
int time_jacobi(timing_result& result, int n) {
    pp::matrix A = random_symmetric_matrix(n);

    jacobi j;
    auto start = std::chrono::high_resolution_clock::now();
    auto [w, V] = j.cyclic(A);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    result.n = n;
    result.time = elapsed.count();
    return 0;
}

// main
int main(int argc, char** argv) {
    // PART 1
    std::cout << "----------PART A----------" << std::endl;
    int n = 3;
    pp::matrix A = random_symmetric_matrix(n);

    std::cout << "-----Checking jacobi-----" << std::endl;
    A.print("A=");
    jacobi j;
    auto [w, V] = j.cyclic(A);
    V.print("V = ");
    w.print("w = ");
    pp::matrix D(w.size(), w.size());
    for (int i=0; i<D.size1(); i++) {
        D[i][i] = w[i];
    }
    D.print("D =");
    pp::matrix VDV = V * D * V.transpose();
    VDV.print("V * D * V^T = A =");

    pp::matrix VVT = V * V.transpose();
    pp::matrix VTV = V.transpose() * V;
    VVT.print("V * V^T = ");
    VTV.print("V^T * V = ");

    // PART 2
    std::cout << "----------PART B----------" << std::endl;
    std::cout << "Build Hamiltonian" << std::endl;
    
    // define default values if inputs are not provided
    double rmax = 10.0;
    double dr   = 0.1;
    std::string wf_file = "";
    std::string diag_file = "";

    // see which inputs are provided and update accordingly
    for(int i = 1; i < argc; i++){
        std::string arg = argv[i];

        if(arg == "-wf" && i + 1 < argc){
            wf_file = argv[++i];
        }
        if(arg == "-diagf" && i + 1 < argc){
            diag_file = argv[++i];
        }
        if(arg == "-rmax" && i + 1 < argc){
            rmax = atof(argv[++i]);
        }
        if(arg == "-dr" && i + 1 < argc){
            dr = atof(argv[++i]);
        }
    }

    // given code to build the hamiltonian
    int npoints = (int)(rmax/dr) - 1;
    pp::vector r(npoints);
    for(int i=0; i<npoints; i++){
        r[i] = dr * (i + 1);
    }
    pp::matrix H(npoints, npoints);
    for(int i=0; i<npoints-1; i++){
        H[i, i]  = -2 * (-0.5 / dr / dr);
        H[i,i + 1]= 1 * (-0.5 / dr / dr);
        H[i + 1, i]= 1 * (-0.5 / dr / dr);
    }
    H[npoints - 1, npoints - 1]= -2 * (-0.5 / dr / dr);
    for(int i=0; i<npoints; i++){
        H[i, i] += -1 / r[i];
    }

    // Diagonalization of H and saving its output to stderr
    std::cout << "Diagonalizing H and finding eigenvalues and eigenvectors" << std::endl;
    auto [wH, VH] = j.cyclic(H);

    std::cerr << rmax << " " << dr << " " << wH[0] << std::endl;

    // saving the wanted wavefunctions in specified file provided by input
    if (wf_file != "") {
        std::ofstream wf(wf_file);
        for(int i = 0; i < npoints; i++){
            wf << r[i];
            for(int k = 0; k < 3; k++){
                float VH_NORM = 1 / (std::sqrt(dr)) * VH[i][k];
                wf << " " << VH_NORM;
            }
            wf << "\n";
        }
        wf.close();
    }

    // PART 3
    std::cout << "----------PART C----------" << std::endl;

    // parrallel processing the matrix diagonalization and saving data into file specified as input
    std::vector sizes = linspace(10, 300, 50);
    int nsizes = sizes.size();
    std::vector<std::thread> threads;
    threads.reserve(nsizes);
    std::vector<timing_result> results(nsizes);

    for (int i=0; i<nsizes; i++) {
        threads.emplace_back(time_jacobi, std::ref(results[i]), sizes[i]);
    }

    for (std::thread& t : threads) {
        t.join();
    }

    if (diag_file != "") {
        std::ofstream df(diag_file);
        for (int i = 0; i < nsizes; i++) {
            df << results[i].n << " " << results[i].time << std::endl;
        }
        df.close();
    }

    return 0;
}