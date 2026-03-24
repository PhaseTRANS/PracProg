#include<iostream>
#include<cmath>
#include<thread>
#include<string>
#include<vector>

struct datum{
    int start, end; 
    double sum;
};

void harm(datum& p) { // harm == harmonic sum
    int start=p.start, end=p.end;
    double sum=0;
    for (int i=start; i<end; i++) {
        sum += 1.0/i;
    }
    p.sum=sum;

}

int main(int argc, char** argv) {
    int nterms=(int)1e9, nthreads=1;
    for(int i=0; i<argc; i++) {
        std::string arg=argv[i];
        if(arg=="-terms" && i+1<argc) {
            nterms=(int)std::stod(argv[++i]);
        }
        if(arg=="-threads" && i+1<argc) {
            nthreads=(int)std::stoi(argv[++i]);
        }
    }
    std::cout <<"terms: "<<nterms<<"\n";
    std::cout <<"threads: "<<nthreads<<"\n";
    
    std::vector<std::thread> threads;  // only variable is created and no constructer
    threads.reserve(nthreads);
    std::vector<datum> data(nthreads);

    for(int i; i<nthreads; ++i) {
        data[i].start=1+(nterms/nthreads)*i;
        data[i].end=1+(nterms/nthreads)*(i+1);
        threads.emplace_back(harm, std::ref(data[i]));
    }
    for(std::thread &thread : threads){
        thread.join();
    }
    double total=0;
    for(datum &d : data){
        total+=d.sum;
    }
    std::cout << "total sum=" << total << std::endl;

}
