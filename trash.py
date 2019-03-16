'''This file contains some methods that were constructed to try gradient descent to find an optimal handle.
I now realize that this was a waste of time, because any handle will do. But I'm keeping this in case I need it for reference. It's probably trash.
def error(alphas, vectors):
    #print alphas
    handle = reduce(lambda x,y: x+y, [a*v for a,v in zip(alphas,vectors)], np.zeros(4))
    #if (np.linalg.norm(handle) - 1) >epsilon:
    #    print "handle not normalized"
    #The next equation appears to be square invariant (handle is a linear function of the alphas). It shouldn't matter.
    error = [((np.dot(handle, vectors[i])/alphas[i])**2*1/(np.dot(handle,handle)) - 4.5)**2 for i in range(len(vectors))]
    error = reduce(lambda x,y: x+y, error, 0)
    return error
def error_derivative(alphas, vectors, k):
    #Returns the partial derivative of the error with respect to alpha_k
    total_sum = 0
    handle = reduce(lambda x,y: x+y, [a*v for a,v in zip(alphas,vectors)], np.zeros(4))
    for i in range(len(vectors)):
        factor1 = 2*((np.dot(handle, vectors[i])/alphas[i])**2/np.dot(handle,handle)-4.5)

        helper_sum1 = 0
        if i != k:
            for j in range(len(vectors)):
                if j!=k:
                    helper_sum1 += alphas[j]/alphas[i]*np.dot(vectors[j],vectors[i])
            #factor2pt1 = 2*np.dot(handle,vectors[i])/alphas[i]*helper_sum1
            factor2pt1 = 2*np.dot(handle,vectors[i])/alphas[i]*np.dot(vectors[k],vectors[i])/alphas[i]*(1/np.dot(handle,handle))
        else:
            factor2pt1 = 0
            for j in range(len(vectors)):
                if j!= k:
                    helper_sum1 += -alphas[j]/alphas[k]**2*np.dot(vectors[j],vectors[k])
            factor2pt1 = 2*np.dot(handle,vectors[i])/alphas[i]*helper_sum1/np.dot(handle,handle)
        helper_sum2 = 0
        for j in range(len(vectors)):
            if j!= k:
                helper_sum2 += alphas[j]*np.dot(vectors[k],vectors[j])
        factor2pt2 = -(np.dot(handle, vectors[i])/alphas[i])**2 * (np.dot(handle,handle)**-2) * (2*alphas[k] + 2*helper_sum2)
        total_sum += factor1*(factor2pt1 + factor2pt2)


    return total_sum
def alphas_from_handle(handle, vectors):
    #Returns a list of alphas so that the handle is their sum
    alphas = [np.dot(handle,v) for v in vectors]
    #print alphas[0]
    #print vectors[0]
    [alphas[i]*vectors[i] for i in range(len(vectors))]
    large_handle = reduce(lambda x,y: x+y, [alphas[i]*vectors[i] for i in range(len(vectors))], np.zeros(4))
    N = np.linalg.norm(large_handle)
    alphas = [a/N for a in alphas]
    return alphas

def renormalize_alphas(alphas):
    #enfoces that the trace must be 1.
    norm = reduce(lambda x,y: x+y, [a**2 for a in alphas], 0)
    normalized_alphas = [a/norm for a in alphas]
    return normalized_alphas
'''
'''
def find_optimal_handle(vectors):
    def adjustments(ws, damping):
        def average_squares(ws):
            inner_prods = [np.dot(handle, w) for w in ws]
            squared_inner_prods = [i**2 for i in inner_prods]
            average_squares = reduce(lambda x,y: x + y, squared_inner_prods, 0) / float(len(ws))
            return average_squares
        average_squares = average_squares(ws)
        adjustments = [(average_squares-np.dot(handle,w))*damping for w in ws]
        return adjustments

    stop = False
    handle = np.array(np.random.rand(4))
    handle = handle/np.linalg.norm(handle)
    ws = map(lambda x: np.dot(handle, x)*x, vectors)
    while(not stop):
        new_ws = map(lambda x: np.dot(handle,x)*x, vectors)
        new_handle = reduce(lambda x,y: x+y, new_ws, np.zeros(4))
        new_handle = new_handle/np.linalg.norm(new_handle)
        thetas = [np.dot(new_handle, v)/np.linalg.norm(w) for v,w in zip(vectors, new_ws)]
        if(all(abs(np.dot(new_handle, vectors[0])/np.linalg.norm(new_ws[0])- np.dot(new_handle, vectors[t])/np.linalg.norm(new_ws[t]))<epsilon for t in range(len(vectors)))):
            stop = True
        else:
            #print np.linalg.norm(new_handle-handle)
            print max(adjustments(ws,0.5))
            handle = new_handle
            ws = [w+x for (x,w) in zip(adjustments(ws,0.5), ws)]
    return new_handle

    '''
#cab_vects = zip(map(lambda x: np.array(x)/np.linalg.norm(x), cab_vects ), map(lambda x: str(x), cab_vects))
