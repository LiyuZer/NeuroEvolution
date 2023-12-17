
class Synapse:
    def __init__(self, id, weight, next_neuron, previous_neuron):
        self.previous_neuron = previous_neuron
        self.next_neuron = next_neuron
        self.weight = weight
        self.activations = [0,0]
        self.id = id
        self.current_time = 0

    def add_activation(self, signal):
        if signal['time'] >= self.current_time: 
            self.activations[0] = (signal['value'])
            

    def activate(self):
        self.current_time = self.current_time + 1
        self.next_neuron.add_activation({'value': self.activations[1] * self.weight, 'time': self.current_time})
        self.activations[1] = self.activations[0]
        self.activations[0] = 0

    
        

    
