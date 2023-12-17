import random
import matplotlib.pyplot as plt
import networkx as nx
from neuron import Neuron
from synapse import Synapse

div = 1000
max_layers = 8
class Serializer:
    def __init__(self) -> None:
        self.nodes = {}
    def add_node(self, node_id, is_input, is_output, attributes):
        self.nodes[node_id/div] = {'connections' : [], 'weights': [], 'val': [is_input, is_output] + attributes}
    def connection(self, node_id, other_id, weight):
        self.nodes[node_id/div]['connections'].append(other_id/div)
        self.nodes[node_id/div]['weights'].append(weight)
    def to_vector(self):
        vector = []
        vector.append(-1)  # Separator
        for node_id, node_data in self.nodes.items():
            vector.append(node_id)
            vector.extend(node_data['val'])  # is_input, is_output, attributes
            for other_id, weight in zip(node_data['connections'], node_data['weights']):
                vector.extend([other_id, weight])
            vector.append(-1)  # Separator
        return vector

class Manager:
    def __init__(self, input_dim=0, hidden_dim_array=[], output_dim=0, array=[]) -> None:
        self.input_neurons = []
        self.output_neurons = []
        self.hidden_neurons = []
        self.synapses = []
        self.id = 1
        self.serializer = Serializer()

        if len(array) == 0:
            # Create neurons
            for i in range(input_dim):
                self.input_neurons.append(Neuron(self.id, False, 0.99, 0.2))
                self.serializer.add_node(self.id, 1, 0, [0.99, 0.3])
                self.id += 1
            for i in range(output_dim):
                self.output_neurons.append(Neuron(self.id, True, 0.99, 0.2))
                self.serializer.add_node(self.id, 0, 1, [0.99, 0.3])
                self.id += 1
            layer_num = 1
            for hidden_layer in hidden_dim_array:
                hidden_layer_neurons = []
                for _ in range(hidden_layer):
                    self.serializer.add_node(self.id, 0, 0, [layer_num/max_layers, 0.99, 0.3])
                    hidden_layer_neurons.append(Neuron(self.id, False, 0.99, 0.2 ))
                    self.id += 1
                self.hidden_neurons.append(hidden_layer_neurons)
                layer_num += 1
            
            # Connect neurons
            self.connect_neurons()
        else:
            nodes = []
            sub_array = []
            for i in array:
                
                if i == -1:
                    nodes.append(sub_array)
                    sub_array = []
                else:
                    sub_array.append(i)
            
            nodes = nodes[1:]
            self.layers = []
            map_node = {}
            for i in nodes:
                size = len(i)
                if i[size - 4] == 1:
                    node = Neuron(self.id, False, i[size - 2], i[size - 1])
                    self.input_neurons.append(node)
                    self.id += 1
                    map_node[self.id] = node
                elif i[size - 3] == 0:
                    node = Neuron(self.id, True, i[size - 2], i[size - 1])
                    self.output_neurons.append(node)
                    self.id += 1
                    map_node[self.id] = node
                else:
                    layer_index = (i[size - 3]) * max_layers
                    # Ensure the layer exists in self.layers
                    while len(self.layers) <= layer_index:
                        self.layers.append([])
                    # Add neuron to the specific layer
                    print(i, "I am here")
                    node = Neuron(self.id, False, i[size - 2], i[size - 1])
                    self.layers[int(layer_index)].append(node)
                    self.id += 1
                    map_node[self.id] = node
            
            for i in nodes:
                weights = i[1:-4]
                size = len(weights)
                for v in range(size):
                    other_node = weights[v] * div
                    strength = weights[v+1]
                    v = v + 1
                    self.create_synapse(map_node[i[0] * div], map_node[other_node], strength)

                




    def connect_neurons(self):
        # Connect input neurons to hidden layers
        for input_neuron in self.input_neurons:
            for layer_index, hidden_layer in enumerate(self.hidden_neurons):
                for hidden_neuron in hidden_layer:
                    if random.random() < 1 / (layer_index + 1):  # Decreasing probability for further layers
                        self.create_synapse_random(input_neuron, hidden_neuron)

        # Recurrent connections in hidden layers and connections from output to hidden
        for layer in self.hidden_neurons:
            for neuron in layer:
                # Recurrent connections within the layer
                for other_neuron in layer:
                    if neuron != other_neuron and random.random() < 0.1:  # Arbitrary probability for recurrent connection
                        self.create_synapse_random(neuron, other_neuron)

                # Connections from output to this neuron in hidden layer
                for output_neuron in self.output_neurons:
                    if random.random() < 0.2:  # Arbitrary probability for connection from output to hidden
                        self.create_synapse_random(output_neuron, neuron)
        

        for output_neuron in self.output_neurons:
            # Select a random neuron from the last hidden layer (or input layer if no hidden layers)
            if self.hidden_neurons:
                connect_from_layer = self.hidden_neurons[-1]
            else:
                connect_from_layer = self.input_neurons

            # Choose a random neuron from the selected layer to connect
            neuron_to_connect = random.choice(connect_from_layer)

            # Create a synapse from the chosen neuron to the current output neuron
            self.create_synapse_random(neuron_to_connect, output_neuron)

        for layer in self.hidden_neurons:
            for neuron in layer:
                if random.random() < 0.1:  # Arbitrary probability for self-recurrent connection
                    self.create_synapse_random(neuron, neuron)


    def create_synapse(self, from_neuron, to_neuron, weight):
        synapse = Synapse(self.id, weight , to_neuron, from_neuron)
        self.synapses.append(synapse)
        from_neuron.forward_synapses.append(synapse)
        # self.serializer.connection(from_neuron.get_id(), to_neuron.get_id(), weight)
        self.id += 1

    def create_synapse_random(self, from_neuron, to_neuron):
        weight = random.uniform(0, 1)
        synapse = Synapse(self.id, weight , to_neuron, from_neuron)
        self.synapses.append(synapse)
        from_neuron.forward_synapses.append(synapse)
        self.serializer.connection(from_neuron.get_id(), to_neuron.get_id(), weight)
        self.id += 1
    def input(self, input):
        count = 0
        for input_neuron in self.input_neurons:
            input_neuron.add_input(input[0])
            count += 1
    def activate(self):
        # Activate all synapses
        for synapse in self.synapses:
            synapse.activate()

        # Activate all neurons
        for neuron_list in [self.input_neurons] + self.hidden_neurons + [self.output_neurons]:
            for neuron in neuron_list:
                neuron.activate()

    def serialize(self):
        return self.serializer.to_vector()
def visualize_network(manager):
    G = nx.DiGraph()

    # Define positions for each layer
    total_layers = 2 + len(manager.hidden_neurons)  # Input, output, and hidden layers
    layer_width = 1.0 / total_layers
    x_positions = {i: (i + 1) * layer_width for i in range(total_layers)}

    # Function to get the value of a neuron
    def get_neuron_value(neuron):
        return neuron.value[1] if hasattr(neuron, 'value') and len(neuron.value) > 1 else "N/A"

    # Add neurons as nodes and determine their positions
    pos = {}
    for i, neuron_list in enumerate([manager.input_neurons] + manager.hidden_neurons + [manager.output_neurons]):
        y_step = 1.0 / (len(neuron_list) + 1)
        for j, neuron in enumerate(neuron_list):
            G.add_node(neuron.id, value=get_neuron_value(neuron))
            pos[neuron.id] = (x_positions[i], (j + 1) * y_step)

    # Add synapses as edges
    for synapse in manager.synapses:
        if synapse.previous_neuron.id == synapse.next_neuron.id:
            # Handle self-recurrent connection
            pos_high = (pos[synapse.previous_neuron.id][0], pos[synapse.previous_neuron.id][1] + 0.1)  # Slightly above the neuron
            G.add_edge(synapse.previous_neuron.id, synapse.next_neuron.id, weight=synapse.weight)
            nx.draw_networkx_edges(G, pos={synapse.previous_neuron.id: pos[synapse.previous_neuron.id], synapse.next_neuron.id: pos_high}, 
                                   edgelist=[(synapse.previous_neuron.id, synapse.next_neuron.id)], 
                                   connectionstyle='arc3,rad=0.2', 
                                   arrows=True)
        else:
            # Normal connection
            G.add_edge(synapse.previous_neuron.id, synapse.next_neuron.id, weight=synapse.weight)

    # Draw the graph
    edge_labels = {(synapse.previous_neuron.id, synapse.next_neuron.id): f"{synapse.weight:.2f}" for synapse in manager.synapses}
    nx.draw(G, pos, with_labels=False, node_size=700, node_color="skyblue", font_size=10, font_weight="bold")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')

    # Display neuron values as node labels
    node_labels = {node: f"{G.nodes[node]['value']:.2f}" for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_color='green')

    plt.show()



# Example usage
array_input = []
for i in range(3):
    array_input.append(random.randint(1,100))

manager = Manager(6, [3,3,3], 6, [])
# manager = Manager(0,0,[],[-1, 0.001, 0.007, 0.7770767823985489, 1, 0, 0.99, 0.3, -1, 0.002, 0.007, 0.03847206874886022, 1, 0, 0.99, 0.3, -1, 0.003, 0.007, 0.7206822514004416, 1, 0, 0.99, 0.3, -1, 0.004, 0.007, 0.7024340755866705, 1, 0, 0.99, 0.3, -1, 0.005, 0.007, 0.3243145406223863, 1, 0, 0.99, 0.3, -1, 0.006, 0, 1, 0.99, 0.3, -1, 0.007, 0.006, 0.4488038075979046, 0, 0, 0.125, 0.99, 0.3, -1])
manager.input(array_input)
print(manager.serialize())
for i in range(10):
    manager.activate()
    visualize_network(manager)

