//
//  Utilities.hpp
//  openCVtest
//
//  Created by Liyu Zerihun on 7/29/22.
//
#include <iostream>
#include <random>
#include <vector>

#ifndef Utilities_hpp
#define Utilities_hpp

extern double learningRate;
using namespace std;

static std::random_device rd;      // obtain a random number from hardware
static std::mt19937 gen(rd());     // seed the generator

static double uniformTest(double x, double y) {

    std::uniform_int_distribution<> distr(x, y);  
    return (double)distr(gen);
}


static double uniformTestRange(double x, double y, int division) {
    float dif=y-x;
    float multiple=dif/division;
    std::uniform_int_distribution<> distr(0, division);  
    return (double)distr(gen)*multiple+x;
}

const int imageWidth = 500;
const int imageHeight = 500;

struct RGB {
    double r;
    double b;
    double g;
    RGB(double red, double blue, double green) {
        r = red;
        b = blue;
        g = green;
    }
    void red(double red) {  // Methods for RGB value modification
        r = red;
    }
    void green(double green) { g = green; }
    void blue(double blue) { b = blue; }
};

// static void seLu(vector<double>* input, vector<double>* output) {
//     for (int i = 0; i < input->size(); i++) {
//         double value = input->at(i);  // the final element of the input matrix
//         if (value < 0) {
//             output->push_back((double)0.001 * value);
//         } else {
//             output->push_back(value);
//         }
//     }
// }

static void sigmoid(vector<double>* input, vector<double>* output) {
    for (int i = 0; i < input->size(); i++) {
        output->push_back(1 / (1 + pow(exp(1), -1 * input->at(0))));
    }
}

static void tanh(vector<double>* input, vector<double>* output) {
    for (int i = 0; i < input->size(); i++) {
        output->push_back((exp(input->at(0)) - exp(-1 * input->at(0))) / (exp(input->at(0)) + exp(-1 * input->at(0))));
    }
}

static void derTanh(double& derivative, double output, double specificOutput, double numberOf, vector<double>& list) {
    derivative = 1 - pow((exp(output) - exp(-1 * output)) / (exp(output) + exp(-1 * output)), 2);
}

static void derReLu(double& derivative, double output, double specificOutput, double numberOf, vector<double>& list) {
    if (output < 0) {
        derivative = 0.001;
    } else if (output == 0) {
        derivative = 0;
    } else {
        derivative = 1;
    }
}

static void reLu(vector<double>* input, vector<double>* output) {
    for (int i = 0; i < input->size(); i++) {
        double value = input->at(i);  // the final element of the input matrix
        if (value <= 0) {
            output->push_back(0.001 * value);
        } else {
            output->push_back(value);
        }
    }
}

static void derSigmoid(double& derivative, double output, double specificOutput, double numberOf,
                       vector<double>& list) {
    derivative = (1 / (1 + pow(exp(1), -1 * output))) * (1 - (1 / (1 + pow(exp(1), -1 * output))));
}

static void derMultiply(double& derivative, double output, double specificOutput, double numberOf,
                        vector<double>& list) {
    double number = 1;
    if (output == 0 || specificOutput == 0) {
        for (int i = 0; i < list.size(); i++) {
            double num = list.at(i);
            if (num != specificOutput) {
                number = number * num;
            }
        }

        derivative = number;
    } else {
        derivative = (double)output / specificOutput;
    }
}

static void derAddition(double& derivative, double output, double specificOutput, double numberOf,
                        vector<double>& list) {
    derivative = 1;
}

static void multiply(vector<double>* input, vector<double>* output) {
    double product = 1;
    for (int i = 0; i < input->size(); i = i + 1) {
        product = product * input->at(i);
    }
    output->push_back(product);
}

static void subtraction(vector<double>* input, vector<double>* output) {
    double difference = 0;
    for (int i = 0; i < input->size(); i = i + 1) {
        difference = difference - input->at(i);
    }
    output->push_back(difference);
}

static void derSubtraction(double& derivative, double output, double specificOutput, double numberOf,
                           vector<double>& list) {
    for (int i = 0; i < list.size(); i = i + 1) {
        if (list.at(i) == specificOutput) {
            if (i % 2 == 0) {
                derivative = 1;
            } else {
                derivative = -1;
            }
        }
    }
}

static void add(vector<double>* list, vector<double>* output) {  // Collective node function
    double sum = 0;

    for (int i = 0; i < list->size(); i = i + 1) {
        sum = list->at(i) + sum;
    }
    output->push_back(sum);
}

static void average(vector<double>* list, vector<double>* output) {  // Collective node function
    double sum = 0;

    for (int i = 0; i < list->size(); i = i + 1) {
        sum = list->at(i) + sum;
    }

    sum = sum / list->size();
    output->push_back(sum);
}

static void not_gate(vector<double>* list, vector<double>* output) {  // Collective node function
    output->push_back(!(list->at(0)));
}

static void derAverage(double& derivative, double output, double specificOutput, double numberOf,
                       vector<double>& list) {  // Collective node function
    derivative = double(1 / numberOf);
}

typedef void (*collectFunc)(
    vector<double>*,
    vector<double>*);  // Collective node, with a vector being accepted as a parameter and a void return statement
typedef void (*collectFuncDerivative)(
    double& derivative, double output, double specificOutput, double number,
    vector<double>& list);  // Collective node, with a vector being accepted as a parameter and a void return statement

static collectFunc collectiveFunc[8] = {&add, &reLu, &multiply, &average, &subtraction, &sigmoid, &tanh, &not_gate};
static collectFuncDerivative derivativeFunc[7] = {&derAddition,    &derReLu,    &derMultiply, &derAverage,
                                                  &derSubtraction, &derSigmoid, &derTanh};


static float tanhyper(float x){
    return (exp(x)-exp(-1*x))/(exp(x)+exp(-1*x));
}
// static void calculateMeanAndVariance(vector<double> vals, double& mean, double& variance) {
//     for (int i = 0; i < vals.size(); i++) {
//         mean = vals.at(i) + mean;
//     }
//     mean = mean / vals.size();

//     for (int i = 0; i < vals.size(); i++) {
//         variance = pow(vals.at(i) - mean, 2) + variance;
//     }

//     variance = variance / vals.size();
// }

#endif /* Utilities_hpp */
