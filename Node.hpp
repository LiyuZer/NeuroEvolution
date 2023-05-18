//
//  Node.hpp
//  openCVtest
//
//  Created by Liyu Zerihun on 8/1/22.
//
#ifndef Node_hpp
#define Node_hpp

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "Network.hpp"
#include "Utilities.hpp"

const double momentum = 0.9;
const int threadCount = 32;
using namespace std;
class Node {
   public:
    std::mutex mut;

    Node(int index, bool is_Value, Network<Node>* n, bool dynamic, int threadID) {
        thread_ID = threadID;
        indexFunction = index;
        is_collective_node = (!is_Value);  // If it is a value then it will switch
        is_value_node = (is_Value);
        is_parameter = false;
        neuralNet = n;
        active = true;
        dynamicallyAllocated = dynamic;
        seen = false;
        derNumber = 0;
        derivativeAtNode = 0;
        //      forwardSeen=false;
        if (dynamic && n != nullptr) {
            n->push(this, threadID);
        }

        previousNumber = 0;
    }

    virtual void activate() = 0;                     // The activation function shall be described in sub classes
    virtual void addInput(double num, Node* prev) {  // The add input function shall be described in the sub classes

        if (prev != nullptr) {
            previous.push_back(prev);
        }
        input.push_back(num);
    }
    virtual void addI(double num) {  // The add input function shall be described in the sub classes
        mut.lock();
        input.push_back(num);

        if (is_value_node ==
            true) {  // if the value node is true then the activated value is the same as teh is_value node
            // I put this here because teh initial creation of an empty node will not have a input and so no activated
            // value for propagation, so once it is a added to the input vector it will be passed to the activated value
            // and so it is ready for propagation
            setActivatedValue(num);
        };
        mut.unlock();
    }

    virtual void addOutput(double num) {  // The add output function shall be described in the sub classes
        mut.lock();
        output.push_back(num);
        activated_value = num;
        mut.unlock();
    }

    virtual void special_propagate() {}
    void addPrev(Node* pointer) {  // The add input function shall be described in the sub classes

        mut.lock();
        previous.push_back(pointer);
        mut.unlock();
    }
    void addNext(Node* pointer) {
        mut.lock();
        next.push_back(pointer);
        mut.unlock();
    }
    vector<Node*>* getNext() { return &next; }
    virtual double getInput(int index) { return input.at(index); }
    vector<double>* getInputVector() { return &input; }
    vector<double>* getIVector() const { return &input; }
    vector<double>* getOVector() const { return &output; }
    vector<double>* getOutputVector() { return &output; }
    void setInputVector(const vector<double> s) { input = s; }
    void setOutputVector(const vector<double> s) { output = s; }

    vector<Node*>* getPrevVector() { return &previous; }
    vector<Node*>* getNextVector() { return &next; }
    int getIndexValue() const { return indexFunction; }
    void setIndex(int s) { indexFunction = s; }
    vector<Node*>* getPVector() const { return &previous; }
    vector<Node*>* getNVector() const { return &next; }
    void setPreviousVector(vector<Node*>& s) { previous = s; }
    void setNextVector(vector<Node*>& s) { next = s; }
    void setDerivativeVector(vector<double> s) { derivative = s; }
    virtual void special_activation(vector<double>& weight, double bias) {}

    void printInfo(Node* p) const {
        if (p->is_parameter == true) {
            cout << "Type: parameter node" << endl;
        } else if (p->is_value_node == true) {
            cout << "Type: value node" << endl;
        }

        else {
            cout << "Type: collection node" << endl;
        }
        cout << "Index Value: " << p->indexFunction << endl;
        cout << "input: " << endl;
        if (p->getIVector()->size() > 0) {
            for (int i = 0; i < p->getIVector()->size(); i = i + 1) {
                cout << p->getInput(i) << endl;
            }
        }
        cout << "output: ";

        for (int i = 0; i < p->getOVector()->size(); i = i + 1) {
            cout << p->getOutput(i) << endl;
        }
        if (p->getOVector()->size() == 0) {
            cout << endl;
        }
        cout << "address: " << p << endl;
        cout << "Previous: " << endl;
        for (int i = 0; i < p->getPVector()->size(); i = i + 1) {
            if (p->getPVector()->at(i) != nullptr) {
                cout << p->getPVector()->at(i) << endl;
            }
        }
        cout << "Derivative: " << endl;
        for (int i = 0; i < p->getDerivative()->size(); i = i + 1) {
            cout << p->getDerivative()->at(i) << endl;
        }

        cout << "Derivative At the node: " << endl;

        cout << p->getDerivativeAtNode() << endl;

        cout << "Next: " << endl;
        for (int i = 0; i < p->getNVector()->size(); i = i + 1) {
            if (p->getNVector()->at(i) != nullptr) {
                cout << p->getNVector()->at(i) << endl;
            }
        }
        cout << endl;
    }
    void printInfoConst(Node* const p) const {
        cout << "value: " << input[0] << endl;
        cout << "address: " << p << endl;
        cout << "Previous: " << endl;
        for (int i = 0; i < p->getPVector()->size(); i = i + 1) {
            if (p->getPVector()->at(i) != nullptr) {
                cout << p->getPVector()->at(i) << endl;
            }
        }

        cout << "Next: " << endl;
        for (int i = 0; i < p->getNVector()->size(); i = i + 1) {
            if (p->getNVector()->at(i) != nullptr) {
                cout << getNVector()->at(i) << endl;
            }
        }
        cout << endl;
    }

    const double getOutput(int index) { return output.at(index); }
    int returnIndex() { return indexFunction; }
    void print(Node* head) {
        printInfo(head);

        if (head->getInputVector()->size() != 0) {
        }

        for (int i = 0; i < head->getNextVector()->size(); i++) {
            if (head->getNextVector()->at(i) != nullptr) {
                print(head->getNextVector()->at(i));
            }
        }
    }

    void printPrev(Node* head) {
        printInfo(head);

        if (head->getInputVector()->size() != 0) {
        }

        for (int i = 0; i < head->getPrevVector()->size(); i++) {
            if (head->getPrevVector()->at(i) != nullptr) {
                printPrev(head->getPrevVector()->at(i));
            }
        }
    }

    void printCount(Node* head, int& x) {
        if (head->is_parameter == true && !head->getSeen()) {
            x++;
            head->setSeen();
        }
        if (head->getInputVector()->size() != 0) {
        }

        for (int i = 0; i < head->getPrevVector()->size(); i++) {
            if (head->getPrevVector()->at(i) != nullptr) {
                printCount(head->getPrevVector()->at(i), x);
            }
        }
    }

    void setActivatedValue(double num) { activated_value = num; }
    double getActivatedValue() { return activated_value; }
    double getAValue() const { return activated_value; }
    void operator>>(Node& right) {
        this->addNext(&right);
        right.addPrev(this);
        if (right.is_collective_node == false) {
            right.addDerivative(1);
        }
    }

    void operator|(Node& right) {
        this->addNext(&right);
        right.addPrev(this);
        right.setInputVector(*this->getIVector());
    }

    void operator>>=(Node& right) {
        this->addNext(&right);
        right.addPrev(this);
        this->propagate();
        if (this->is_collective_node && right.is_value_node) {
            right.getDerivative()->push_back(1);
        }
    }
    virtual void propagate() {
        vector<Node*>* nextVector = this->getNextVector();
        for (int i = 0; i < nextVector->size(); i = i + 1) {
            if (nextVector->at(i) != nullptr) {
                if (is_value_node || is_parameter) {
                    // This is here because during activation for collective node the value is already updated (output
                    // value), thus it only takes place if this node is a value node
                    // Whenever a computation is made an output is added to the vector, note that the output vector is
                    // very important for final calculations as it corresponds to the derivative vector
                    this->addOutput(this->activated_value);
                }
                nextVector->at(i)->addI(activated_value);
            }
        }
    }
    bool get_is_collective_node() const { return is_collective_node; }
    bool get_is_parameter() const { return is_parameter; }
    bool get_is_value_node() const { return is_value_node; }
    void setTrueParameter() { is_parameter = true; }
    void setActiveFalse() const { active = false; }
    bool getActive() { return active; }
    Network<Node>* getNetworkPointer() const { return neuralNet; }
    bool isDynamic() const { return dynamicallyAllocated; }
    void dynamicMemory() { dynamicallyAllocated = true; }
    void setSeen() { seen = true; }
    bool getSeen() { return seen; }

    void addDerNumber() { derNumber = derNumber + 1; }
    int getDerNumber() { return derNumber; }
    void zeroDerNumber() { derNumber = 0; }

    void aggregate(double val) {
        if (is_parameter || is_value_node) {
            getInputVector()->at(0) = getInputVector()->at(0) + val;
        } else {
            getInputVector()->push_back(val);
        }
    }
    virtual ~Node() {}

    void updateNode(Node* p, double learningParam) {
        if (p->is_parameter == true) {
            p->getInputVector()->at(0) = p->getInputVector()->at(0) - (double)p->getDerivativeAtNode() * (learningRate);
        }
        for (int i = 0; i < p->getPrevVector()->size(); i = i + 1) {
            updateNode(p->getPVector()->at(i), learningParam);
        }
    }

    void backProp(Node* p, double der, Node* target, int& count, vector<Node*>& endNodes) {
        if (p->is_parameter == true) {
            if (p->getSeen()) {
            } else {
                count++;
                p->setSeen();
                endNodes.push_back(p);
            }
            p->setDerivativeAtNode((double)p->getDerivativeAtNode() + der);
        }

        else {
            p->addDerNumber();

            if (p->getSeen()) {  // if it hasn't been and has not prev vector add

            } else if (p->getPrevVector()->size() == 0) {
                endNodes.push_back(p);
                p->setSeen();
            }

            if (p->getDerNumber() >= p->getNextVector()->size()) {
                double val = 0;

                if (p->getNextVector()->size() == 0) {
                    val = 1;
                } else {
                    if (p->getDerNumber() == p->getNextVector()->size()) {
                        p->setDerivativeAtNode(der + p->getDerivativeAtNode());
                        val = p->getDerivativeAtNode();
                    }
                }
                for (int i = 0; i < p->getPrevVector()->size(); i = i + 1) {
                    double backPropval = val * p->getDerivative()->at(i);
                    backProp(p->getPrevVector()->at(i), backPropval, target, count, endNodes);
                }

            } else {
                p->setDerivativeAtNode(der + p->getDerivativeAtNode());
            }
        }
    }

    int getPreviousNumber() { return previousNumber; }
    void setPrevNumber(int num) { previousNumber = num; }

    void setDerivativeAtNode(double d) { derivativeAtNode = d; }
    double getDerivativeAtNode() { return derivativeAtNode; }

    void setThreadID(int d) const { thread_ID = d; }
    double getThreadID() const { return thread_ID; }

    virtual void specialActivation() = 0;
    vector<double>* getDerivative() const { return &derivative; }

    void addDerivative(double d) const { derivative.push_back(d); }
    string getHash() const { return hash; }
    void setHash(string s) { hash = s; }

   private:
    mutable vector<double> input;  // input of the node
    double activated_value;
    mutable vector<double> output;  // outputted function
    int indexFunction;               // Index function for array of functions in utility
    mutable vector<Node*> previous;
    mutable vector<Node*> next;
    bool is_value_node;
    bool is_collective_node;
    bool is_parameter;
    mutable vector<double> derivative;
    Network<Node>* neuralNet;
    mutable bool active;
    bool dynamicallyAllocated;
    bool seen;  // Test variable;
    int derNumber = 0;
    mutable vector<double> previousDerivatives;
    string specialName;
    double derivativeAtNode;
    mutable int thread_ID;
    //   bool forwardSeen;
    //    int forwardNumber;
    int previousNumber;
    string hash;
};

class collection_node : public Node {
   public:
    collection_node(int index, Network<Node>* n, bool dynamic, int threadID)
        : Node(index, false, n, dynamic, threadID) {}

    collection_node(Network<Node>* n, bool dynamic, int threadID) : Node(0, false, n, dynamic, threadID) {}

    collection_node(collection_node const& right)
        : Node(0, false, right.getNetworkPointer(), right.isDynamic(), right.getThreadID()) {
        right.setActiveFalse();

        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setIndex(right.getIndexValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setThreadID(right.getThreadID());
        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }
    }

    collection_node& operator=(collection_node& right) {
        mut.lock();
        right.setActiveFalse();
        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setIndex(right.getIndexValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setThreadID(right.getThreadID());

        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }
        mut.unlock();

        return *this;
    }

    void activate() override {
        collectiveFunc[getIndexValue()](getInputVector(), getOutputVector());
        setActivatedValue(getOutput(0));
        propagate();
        for (int i = 0; i < getPrevVector()->size(); i++) {
            double temp;
            derivativeFunc[getIndexValue()](temp, getOutput(0), getPrevVector()->at(i)->getOutput(0),
                                            getInputVector()->size(), *getIVector());
            getDerivative()->push_back(temp);
        }

        getInputVector()->clear();
        getInputVector()->push_back(
            1);  // Once there is an activation the input doesn't matter, take this out for regular use
    }

    void special_propagate() override {
        vector<Node*>* nextVector = this->getNextVector();
        for (int i = 0; i < nextVector->size(); i = i + 1) {
            if (nextVector->at(i) != nullptr) {
                this->addOutput(this->getInputVector()->at(0));  // Whenever a computation is made
            }

            if (nextVector->at(i)->get_is_collective_node()) {
                nextVector->at(i)->addI(getOutputVector()->at(0));
            } else if (nextVector->at(i)->get_is_parameter()) {
            } else {
                nextVector->at(i)->getInputVector()->at(0) =
                    nextVector->at(i)->getInputVector()->at(0) + (getOutputVector()->at(0));
            }
        }
    }
    void special_activation(vector<double>& weight , double bias) {
        // collectiveFunc[getIndexValue()](getInputVector(), getOutputVector());
        // setActivatedValue(getOutput(0));
        // int temp = getInputVector()->size();
        // special_propagate();
        // getInputVector()->erase(getInputVector()->begin(), getInputVector()->begin() + temp);
        // getInputVector()->push_back(1);
        // getOutputVector()->clear();
    }

    void specialActivation() override {};
};

class value_node : public Node {
   public:
    value_node(double value, Node* prev, Network<Node>* n, bool dynamic, int threadID)
        : Node(0, true, n, dynamic, threadID) {
        addInput(value, prev);
        setActivatedValue(value);
    }
    value_node(bool val, Network<Node>* n, bool dynamic, int threadID) : Node(0, true, n, dynamic, threadID) {}

    value_node(Network<Node>* n, bool dynamic, int threadID) : Node(0, true, n, dynamic, threadID) {}

    value_node(const value_node& right)
        : Node(0, true, right.getNetworkPointer(), right.isDynamic(), right.getThreadID()) {
        right.setActiveFalse();
        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setThreadID(right.getThreadID());

        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }
    }

    value_node& operator=(const value_node& right) {
        right.setActiveFalse();
        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setThreadID(right.getThreadID());

        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }

        return *this;
    }

    virtual value_node& operator+(value_node& right) {
        value_node* sum_pointer = new value_node(this->getInput(0) + right.getInput(0), &right,
                                                 right.getNetworkPointer(), true, this->getThreadID());
        sum_pointer->addPrev(this);

        this->addNext(sum_pointer);
        right.addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        right.addOutput(right.getInput(0));
        sum_pointer->getDerivative()->push_back(1);
        sum_pointer->getDerivative()->push_back(1);

        return *sum_pointer;
    }

    value_node& operator&(double right) {
        value_node* sum_pointer =
            new value_node(pow(exp(1), this->getInput(0)), this, this->getNetworkPointer(), true, this->getThreadID());
        this->addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        sum_pointer->getDerivative()->push_back(pow(exp(1), this->getInput(0)));
        return *sum_pointer;
    }

    value_node& operator%(double right) {
        value_node* sum_pointer =
            new value_node(log(this->getInput(0)), this, this->getNetworkPointer(), true, this->getThreadID());

        this->addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        sum_pointer->getDerivative()->push_back((double)1 / this->getInput(0));
        return *sum_pointer;
    }

    value_node& operator-(value_node& right) {
        value_node* sum_pointer = new value_node(this->getInput(0) - right.getInput(0), &right,
                                                 right.getNetworkPointer(), true, this->getThreadID());
        sum_pointer->addPrev(this);
        this->addNext(sum_pointer);
        right.addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        right.addOutput(right.getInput(0));
        sum_pointer->getDerivative()->push_back(-1);
        sum_pointer->getDerivative()->push_back(1);

        return *sum_pointer;
    }

    value_node& operator/(value_node& right) {
        value_node* sum_pointer = new value_node(this->getInput(0) / right.getInput(0), &right,
                                                 right.getNetworkPointer(), true, this->getThreadID());
        sum_pointer->addPrev(this);

        this->addNext(sum_pointer);
        right.addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        right.addOutput(right.getInput(0));
        sum_pointer->getDerivative()->push_back(-1 * this->getInput(0) * ((double)1 / pow(right.getInput(0), 2)));
        sum_pointer->getDerivative()->push_back((double)1 / right.getInput(0));

        return *sum_pointer;
    }
    virtual value_node& operator*(value_node& right) {
        value_node* sum_pointer = new value_node(this->getInput(0) * right.getInput(0), this, right.getNetworkPointer(),
                                                 true, this->getThreadID());
        sum_pointer->addPrev(&right);
        this->addNext(sum_pointer);
        right.addNext(sum_pointer);

        this->addOutput(this->getInput(0));
        right.addOutput(right.getInput(0));
        sum_pointer->addDerivative(right.getInput(0));
        sum_pointer->addDerivative(this->getInput(0));

        return *sum_pointer;
    }

    value_node& operator^(double right) {
        value_node* sum_pointer = new value_node((double)pow(this->getInput(0), right), this, this->getNetworkPointer(),
                                                 true, this->getThreadID());

        this->addNext(sum_pointer);
        this->addOutput(this->getInput(0));
        sum_pointer->getDerivative()->push_back((double)(right)*pow(this->getInput(0), right - 1));
        return *sum_pointer;
    }
    virtual double getInput(int index) { return Node::getInput(index); }

    void activate() {
        setActivatedValue(getInput(0));
        propagate();
        getInputVector()->clear();
        getInputVector()->push_back(0);
    }

    void special_activation(vector<double>& weight, double bias) {
        for(int i=0; i<weight.size(); i++){
        getOutputVector()->push_back(rellu(getInputVector()->at(0)*weight.at(i)+bias));
        }
        double temp = getInputVector()->at(0);
        special_propagate();
        getInputVector()->at(0) = getInputVector()->at(0) - temp;
        getOutputVector()->clear();
    }
    void specialActivation() { setActivatedValue(getInput(0)); }
    void special_propagate() {
        
        vector<Node*>* nextVector = this->getNextVector();
        for (int i = 0; i < nextVector->size(); i = i + 1) {


            if (nextVector->at(i)->get_is_collective_node()) {
                nextVector->at(i)->addI(getInputVector()->at(i));
            } else if (nextVector->at(i)->get_is_parameter()) {
            } else {
                nextVector->at(i)->getInputVector()->at(0) =
                    nextVector->at(i)->getInputVector()->at(0) + getOutputVector()->at(i);
            }
        }
    }
};

class parameter_node : public value_node {
   public:
    parameter_node(double value, Node* prev, Network<Node>* n, bool dynamic, std::shared_ptr<double[]> pre,
                   int threadID)
        : value_node(n, dynamic, threadID) {
        if (pre != nullptr) {
            previousDerivatives = pre;
            pre.reset();
            param = &previousDerivatives[threadCount];
            increment = &previousDerivatives[threadCount + 1];
            addInput(*param, prev);
        } else {
            std::shared_ptr<double[]> shared(new double[threadCount + 2]());
            previousDerivatives = shared;
            shared.reset();
            for (int i = 0; i < threadCount; i++) {
                previousDerivatives[i] = 0;
            }
            param = &previousDerivatives[threadCount];
            *param = value;
            increment = &previousDerivatives[threadCount + 1];
        }
        setActivatedValue(*param);
        setTrueParameter();
        addInput(*param, prev);
    }
    parameter_node(Network<Node>* n, bool dynamic, shared_ptr<double[]> pre, int threadID)
        : value_node(true, n, dynamic, threadID) {
        setTrueParameter();
        if (pre != nullptr) {
            previousDerivatives = pre;
            pre.reset();
            param = &previousDerivatives[threadCount];
            increment = &previousDerivatives[threadCount + 1];
            addInput(*param, nullptr);

        } else {
            std::shared_ptr<double[]> shared(new double[threadCount + 2]());
            previousDerivatives = shared;
            shared.reset();
            for (int i = 0; i < 6; i++) {
                previousDerivatives[i] = 0;
            }
            param = &previousDerivatives[threadCount];
            increment = &previousDerivatives[threadCount + 1];
            addInput(*param, nullptr);
        }
    }

    parameter_node(Network<Node>* n, bool dynamic, int threadID) : value_node(true, n, dynamic, threadID) {
        setTrueParameter();
        previousDerivatives = nullptr;
    }

    parameter_node(double value, Network<Node>* n, bool dynamic, int threadID)
        : value_node(true, n, dynamic, threadID) {
        setTrueParameter();
        previousDerivatives = nullptr;
        addI(value);
    }

    parameter_node(const parameter_node& right)
        : value_node(right.getNetworkPointer(), right.isDynamic(), right.getThreadID()) {
        right.setActiveFalse();
        setTrueParameter();
        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setPrevDerivative(right.getPrevDerivative());
        this->setThreadID(right.getThreadID());
        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }
        //            cout<<" Being converted to: "<<endl;
        //            this->printInfo(this);
    }

    //    value_node operator=(value_node& right){
    //        this->addNext(&right);
    //        right.addPrev(this);
    //
    //        return *this;
    // }

    parameter_node& operator=(parameter_node& right) {
        setTrueParameter();
        right.setActiveFalse();
        this->setInputVector(*right.getIVector());
        this->setOutputVector(*right.getOVector());
        this->setNextVector(*right.getNVector());
        this->setPreviousVector(*right.getPVector());
        this->setActivatedValue(right.getAValue());
        this->setDerivativeVector(*right.getDerivative());
        this->setPrevDerivative(right.getPrevDerivative());
        this->setThreadID(right.getThreadID());

        for (int i = 0; i < right.getPVector()->size(); i++) {
            if (right.getPVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getPVector()->at(i)->getNextVector()->size(); m = m + 1) {
                    if (right.getPVector()->at(i)->getNVector()->at(m) == &right) {
                        right.getPVector()->at(i)->getNVector()->at(m) = this;
                    }
                }
            }
        }

        for (int i = 0; i < right.getNVector()->size(); i++) {
            if (right.getNVector()->at(i) != nullptr) {
                for (int m = 0; m < right.getNVector()->at(i)->getPVector()->size(); m = m + 1) {
                    if (right.getNVector()->at(i)->getPrevVector()->at(m) == &right) {
                        right.getNVector()->at(i)->getPrevVector()->at(m) = this;
                    }
                }
            }
        }
        //            cout<<" Being converted to: "<<endl;
        //            this->printInfo(this);
        return *this;
    }

    void setParamValue(double val) {
        *param = val;
        this->getIVector()->at(0) = *param;
    }
    void activate() override {
        setActivatedValue(*param);
        propagate();
    }

    void special_activation(vector<double>& weight, double bias)  {
        // setActivatedValue(getInputVector()->at(0));
        // special_propagate();
    }

    virtual double getInput(int index) override { return *param; }

    void special_propagate() override {
        vector<Node*>* nextVector = this->getNextVector();
        for (int i = 0; i < nextVector->size(); i = i + 1) {
            if (nextVector->at(i) != nullptr) {
                this->addOutput(this->getInputVector()->at(0));  // Whenever a computation is made
            }

            if (nextVector->at(i)->get_is_collective_node()) {
                nextVector->at(i)->addI(getInputVector()->at(0));
            } else if (nextVector->at(i)->get_is_parameter()) {
            } else {
                nextVector->at(i)->getInputVector()->at(0) =
                    nextVector->at(i)->getInputVector()->at(0) + (getInputVector()->at(0));
            }
        }
    }

    void propagate() override {
        vector<Node*>* nextVector = this->getNextVector();
        for (int i = 0; i < nextVector->size(); i = i + 1) {
            if (nextVector->at(i) != nullptr) {
                // This is here because during activation for collective node the value is already updated(output
                // value), thus it only takes place if this node is a value node
                this->addOutput(*param);  // Whenever a computation is made an output is added to the vector, note
                // that the output vector is very important for final calculations as it corresponds to the derivative
                // vector
                nextVector->at(i)->addI(*param);
            }
        }
    }

    void update(bool time, int threadID) {
        if (time) {
            double sum = 0;
            for (int i = 0; i < threadCount; i++) {
                sum = this->previousDerivatives[i] + sum;
            }
            double average = sum / threadCount;
            double newWeight = *param - ((double)average * (learningRate) + momentum * (*increment));
            *param = newWeight;
            *increment = (average * (learningRate) + momentum * (*increment));

            for (int i = 0; i < threadCount; i++) {
                previousDerivatives[i] = 0;
            }
        } else {
            this->previousDerivatives[threadID] = ((double)this->getDerivativeAtNode());
        }
    }

    shared_ptr<double[]> getPrevDerivative() const { return previousDerivatives; }

    void setPrevDerivative(shared_ptr<double[]> pre) {
        if (pre == nullptr) {
        } else {
            previousDerivatives = pre;
            pre.reset();
            param = &previousDerivatives[threadCount];
            increment = &previousDerivatives[threadCount + 1];
        }
    }

   private:
    mutable shared_ptr<double[]> previousDerivatives;
    mutable double* param;
    mutable double* increment;
};

#endif /* Node_hpp */
