import itertools

def compatible_list(vertex, list_of_verticies):
    #Finds the neighbors of vertex in list_of_verticies.
    a = len(vertex)/2 #Same a as above
    neighbors = []
    for i in range(a):
        neighbors += [v for v in list_of_verticies if v[i] != vertex[i] and v[i+a]==vertex[i+a]]
        return neighbors

if __name__=='__main__':
    '''
    There are a qubits distributed among a people, each of which can perform b binary experiments.
    When a=b=2, we have Bell's experiment.
    This function returns the exclusivity graph of the experiment.
    The verticies are experiments and their outcomes represented as a string of 1's & 0's, (for the outcomes) followed by a string of length a with letters [0..b]
    They are connected if the two outcomes could not have occured together.
    More precisely, they are connected if the two verticies share an experiment which yeilds different measurements.
    '''
    a = 2
    b = 3
    experiment_choices = range(b)
    experiment_choices = [str(exp) for exp in experiment_choices] #convert to string
    total_experiments = itertools.product(experiment_choices, repeat=a)
    total_outcomes = itertools.product(['0','1'], repeat=a)
    verticies = [''.join(o+e) for o,e in itertools.product(total_outcomes, total_experiments)]
    graph = {node: compatible_list(node, verticies) for node in verticies}
    print graph
