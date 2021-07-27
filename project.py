import math
import matplotlib.pyplot as plt
import numpy as np


########----------1-------------

def minimum_squares_method(function,base,sectionBegin,sectionEnds):
    numOrgansBase=len(base)
    #init the matrix to be with zeros at len n*n
    matpairs=np.zeros((numOrgansBase,numOrgansBase))
    vectorResult=[]
    #clculate the vectorResult
    i=0
    while i<numOrgansBase:
        vectorResult.append(numerical_integral_calculation((lambda x: function(x)*base[i](x)),sectionBegin,sectionEnds))
        i+=1
    print(vectorResult)
    #calculate the matrix
    for i in range(numOrgansBase):
        for j in range(numOrgansBase):
            matpairs[i][j]=numerical_integral_calculation((lambda x: base[i](x)*base[j](x)),sectionBegin,sectionEnds)
    print(matpairs)
    vectorAlfas=calsystemGauss(matpairs,vectorResult)

    ###the new func
    fixfunction = lambda x: vectorAlfas[0] * base[0](x) + vectorAlfas[1] * base[1](x) + vectorAlfas[2] * base[2](x)

    error=calcError(function,fixfunction,sectionBegin,sectionEnds)

    plot_graphs(function, fixfunction)

    return vectorAlfas,error





def calcError(realFunction, fixFunction, a, b):
    g = lambda x: realFunction(x) - fixFunction(x)
    g1 = lambda x: g(x) * g(x)
    error = numerical_integral_calculation(g1, a, b, 50,1)
    return math.sqrt(error)




########---------------2-------------

def numerical_integral_calculation(f,sectionBegin,sectionEnds,n=100,method=2):
    res=0
    if method==1:
        res=trapez(f,sectionBegin,sectionEnds,n)
    else:
        res=simpson(f,sectionBegin,sectionEnds,n)
    return res



def trapez(f,a,b,n):
    h=(b-a)/float(n)
    sigma=0
    for i in range(1,n):
        xi = a + (i * h)
        sigma+=f(xi)
    res=(h/2)*(f(a)+(2*sigma)+f(b))
    return res


def simpson(f,a,b,n):
    sigmaEvenPoints=0.0
    sigmaOddPoints=0.0
    h = (b - a) / n
    for i in range(1,int((n/2)-1)):
        xtwoi = a + (2*i * h)
        sigmaEvenPoints+=f(xtwoi)
    for i in range(1,int((n/2))):
        xtwoi = a + (2 * i * h)
        sigmaOddPoints+=f(xtwoi-1)
    res=(h/3)*(f(a)+(2*sigmaEvenPoints)+(4*sigmaOddPoints)+f(b))
    return res



##########--------------3--------------

def calsystemGauss(matpairs,vectorResult):
    vectorAlfas=np.zeros(len(vectorResult))
    iteration=25
    n=len(vectorResult)
    for k in range(iteration):
        for i in range(n):
            for j in range(n):
                if(i != j):
                    vectorAlfas[i]=(vectorResult[i]-(matpairs[i][j]*vectorAlfas[j]))/matpairs[i][i]
    return vectorAlfas



def graph(func, x_range, cl='r--'):
    y_range=[]
    for x in x_range:
        y_range.append(func(x))
    plt.plot(x_range, y_range, cl)
    return

def plot_graphs(f, ff):
    rs=1.0
    r=np.linspace(-rs*np.pi,rs*np.pi,80)
    graph(ff,r,cl='r-')
    graph(f,r,cl='g--')
    plt.axis('equal')
    plt.show()





if __name__ == '__main__':
    a = -math.pi
    b = math.pi

    # Fourier serias
    a1 = lambda x: 1 / math.sqrt(2 * math.pi)
    a2 = lambda x: (1 / math.sqrt(math.pi)) * math.cos(x)
    a3 = lambda x: (1 / math.sqrt(math.pi)) * math.cos(x)
    base1=[a1,a2,a3]
    f = lambda t: abs(t)
    alfas,error=minimum_squares_method(f,base1,a,b)
    print(f'\nERROR: {error}')
    print(alfas)

    # a0 = lambda x: 1
    # a1 = lambda x: x
    # a2 = lambda x: x ** 2
    # a3 = lambda x: x ** 4
    # base1 = [a0, a1, a2, a3]
    # f = lambda t: abs(t)
    # alfas, error = minimum_squares_method(f, base1, a, b)
    # print(f'\nERROR: {error}')
    # print(alfas)

    a1 = lambda x: 1
    a2 = lambda x: math.cos(x)
    a3 = lambda x: math.sin(x)
    base1 = [a1, a2, a3]
    f = lambda t: t*t
    alfas, error = minimum_squares_method(f, base1, a, b)
    print(f'\nERROR: {error}')
    print(alfas)


