#include<iostream>
#include<cmath>
#include<vector>
#include<functional>
#include<fstream>
#include"matrix.hpp"
#include"newton.hpp"


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

struct ann{
    int n; /* number of hidden neurons */
    std::function<double(double)> f = [](double x) { /* activation function */
        return x*std::exp(-x*x);
    };

    pp::vector p; /* network parameters: 3 per hidden neuron (center, width, weight) */
    ann(int n) : n(n) {
        p.resize(3*n);
        for(int i=0; i<n; i++) {
            p[3*i] = -1.0 + 2.0*i/(n-1);    // center
            p[3*i+1] = 1.0;                 // width
            p[3*i+2] = 1.0;                 // weight
        }
    };

    // derivatives and anit derivative of activation 
    std::function<double(double)> df  = [](double x){ return std::exp(-x*x) * (1 - 2*x*x); };
    std::function<double(double)> ddf = [](double x){ return std::exp(-x*x) * x * (4*x*x - 6); };
    std::function<double(double)> F   = [](double x){ return -0.5 * std::exp(-x*x); };

    double forward(double x) const {
        double y = 0;
        for (int i=0; i<n; i++) {
            double center = p[3*i], width = p[3*i+1], weight = p[3*i+2];
            double z = (x - center) / width;
            y += weight * f(z);
        }
        return y;
    }

    // first derivative
    double derivative(double x) const {
        double y = 0;
        for (int i=0; i<n; i++) {
            double center = p[3*i], width = p[3*i+1], weight = p[3*i+2];
            double z = (x - center) / width;
            y += weight * df(z) * (1.0 / width);
        }
        return y;
    }

    // second derivative
    double second_derivative(double x) const {
        double y = 0;
        for (int i=0; i<n; i++) {
            double center = p[3*i], width = p[3*i+1], weight = p[3*i+2];
            double z = (x - center) / width;
            y += weight * ddf(z) * (1.0 / (width*width));
        }
        return y;
    }

    // antiderivative
    double antiderivative(double x, double x0=0.0) const {
        double y = 0;
        for (int i=0; i<n; i++) {
            double center = p[3*i], width = p[3*i+1], weight = p[3*i+2];
            double z  = (x  - center) / width;
            double z0 = (x0 - center) / width;
            y += weight * width * (F(z) - F(z0));
        }
        return y;
    }

    int train_num(pp::vector xs, pp::vector ys, int epochs=100){
        /* train the network to interpolate the given table {x,y} */
        // function for numerical training does not work good for netwrok training
        std::function<double(pp::vector)> cost_function = [this, &xs, &ys](pp::vector params) {
            double cost = 0;
            for (int i=0; i<(int)xs.size(); i++) {
                double y_pred = 0;
                for (int j=0; j<n; j++) {
                    double center = params[3*j];
                    double width = params[3*j+1];
                    double weight = params[3*j+2];
                    y_pred += f((xs[i] - center) / std::abs(width)) * weight;
                }
                double error = y_pred - ys[i];
                cost += error * error;
            }
            return cost;
        };
        p = pp::newton(cost_function, p, 0.001, epochs);
        return 0;
    };

    int train(pp::vector xs, pp::vector ys, int epochs=100, double lr=0.01) {
        /* train the network to interpolate the given table {x,y} */
        // training function that usese gradient descent and backpropagation instead
        // since numerical gradient estimation doesn't work
        
        // derivative of f(x) = x*exp(-x^2) is f'(x) = exp(-x^2)*(1 - 2*x^2)
        auto df = [](double x) {
            return std::exp(-x*x) * (1 - 2*x*x);
        };

        for (int epoch=0; epoch<epochs; epoch++) {
            pp::vector grad(p.size());
            
            // for each epoch do:
            // -evaulate the network prediction
            // -calculate the cost
            // -use backpropagation with analytical graidents
            // -make the optimization step
            for (int i=0; i<(int)xs.size(); i++) {
                // forward pass
                double y_pred = 0;
                for (int j=0; j<n; j++) {
                    double center = p[3*j];
                    double width  = p[3*j+1];
                    double weight = p[3*j+2];
                    y_pred += f((xs[i] - center) / width) * weight;
                }

                double error = y_pred - ys[i];  // residual

                // backward pass
                for (int j=0; j<n; j++) {
                    double center = p[3*j];
                    double width  = p[3*j+1];
                    double weight = p[3*j+2];
                    double z      = (xs[i] - center) / width;
                    double fz     = f(z);
                    double dfz    = df(z);

                    // analytic gradients of graph
                    grad[3*j]   += 2 * error * weight * dfz * (-1.0 / width);  
                    grad[3*j+1] += 2 * error * weight * dfz * (-z / width);    
                    grad[3*j+2] += 2 * error * fz;                              
                }
            }

            // gradient descent step
            for (int k=0; k<(int)p.size(); k++) {
                p[k] -= lr * grad[k];
            }
        }
        return 0;
    };
};

int main() {
    // PART A
    std::cout << "----------PART A-----------" << std::endl;
    auto g = [](double x) {
        return std::cos(5 * x - 1) * std::exp(- x * x);
    };
    // sample 50 points
    std::vector xs_std = linspace(-1, 1, 50);
    pp::vector xs(xs_std.size());
    pp::vector ys(xs.size());
    for (int i=0; i<(int)xs.size(); i++) {
        xs[i] = xs_std[i];
        ys[i] = g(xs[i]);
    }

    // initialize network and train
    ann network(10);
    network.train(xs, ys, 1000);

    // save the result to file for plotting
    std::vector xs_plot = linspace(-1, 1, 200);
    std::ofstream network_out("network.dat");
    for (int i=0; i<(int)xs_plot.size(); i++) {
        if (i<(int)xs.size()) {
            network_out << xs_plot[i] << " " 
            << network.forward(xs_plot[i]) << " "
            << network.derivative(xs_plot[i]) << " "            // PART B
            << network.second_derivative(xs_plot[i]) << " "     // PART B
            << network.antiderivative(xs_plot[i]) << " "        // PART B
            << g(xs_plot[i]) << " "
            << xs[i] << " "
            << ys[i] << std::endl; 
        }
        else {
            network_out << xs_plot[i] << " " 
            << network.forward(xs_plot[i]) << " "
            << network.derivative(xs_plot[i]) << " "            // PART B
            << network.second_derivative(xs_plot[i]) << " "     // PART B
            << network.antiderivative(xs_plot[i]) << " "        // PART B
            << g(xs_plot[i]) << std::endl;
        }
    }
    network_out.close();

    // PART B
    std::cout << "----------PART B-----------" << std::endl;
    // PART B is implememted in the class above and printed in the same file as PART A

    // PART C
    std::cout << "----------PART C-----------" << std::endl;
    

    return 0;
}