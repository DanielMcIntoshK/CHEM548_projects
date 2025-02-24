import matplotlib.pyplot as plt
import numpy as np

def graph(filename, n,ymin,ymax,shouldshow=True):
    f=open(filename)
    
    lines=f.readlines()

    vs=lines[0].split()
    evals=lines[1].split()

    x=np.linspace(-5,5,100);

    V=[]
    for i in range(100):
        V.append(float(vs[i]))
    plt.plot(x,V)

    for s in range(n):
        y=[]
        for i in range(100):
            y.append(float(lines[i+2].split()[s])+float(evals[s]))

        plt.plot(x,y)
        plt.hlines(float(evals[s]),-5,5,linestyles='dashed')
    plt.ylim(ymin,ymax)
    if(shouldshow):
        plt.show()

def main():
    graph("../Output/Pib", 4,-.5,2)
    graph("../Output/Fin", 4,-.5,3)
    graph("../Output/Pbrec", 4,-.5,2)
    graph("../Output/Harm", 4,-.5,5)
    graph("../Output/Morse1",4,-.5,10,False)
    graph("../Output/Morse2",4,-.5,10)

main()

