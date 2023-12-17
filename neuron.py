
from typing import Any
import numpy as np

class Neuron:
    def __init__(self, id, is_output,decay_rate, threshold ):
        self.id = id
        self.forward_synapses = []
        self.value = [0,0]
        self.decay_rate = decay_rate
        self.threshold = threshold
        self.current_time = 0
        self.is_output = is_output


    def add_input(self, signal):
        self.value[0] = signal

    def get_id(self):
        return self.id
    
    def add_activation(self, signal):
        if self.current_time > signal['time']:  
            self.value[1] += signal['value']
        else:
            self.value[0] += signal['value']

    import numpy as np

    def should_fire(self, lam):
        """
        Determines whether to fire (1) or not (0) based on a Poisson distribution.

        Parameters:
        lam (float): The lambda (average rate of occurrence) for the Poisson distribution.

        Returns:
        int: 1 if fire, 0 otherwise.
        """
        # Generate a Poisson-distributed value
        value = np.random.poisson(lam)

        # Determine to fire or not based on the value
        return 1 if value > 0 else 0


    def activate(self):
        self.current_time +=1
        if self.value[1] > self.threshold and self.should_fire(0.4):
            activated_value = max(0.1 * self.value[1], self.value[1])
            for synapse in self.forward_synapses:
                synapse.add_activation({'value': activated_value, 'time': self.current_time})
            if not self.is_output:
                self.value[1] = 0
        else:
            self.value[1] *= self.decay_rate

        self.value[1] = self.value[1] + self.value[0]
        self.value[0] = 0

        
        


## Just a rant here, so for now each object in the neural circuit will have a current time that must be synced up with the current time. Thus each 
## Neuron, synapse and value must be synced up to function, the great thing about this is that parts of the circuit can work independently with other parts
## As they can follow the absolute time 