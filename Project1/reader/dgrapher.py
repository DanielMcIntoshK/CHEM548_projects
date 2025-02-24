import matplotlib.pyplot as plt
import numpy as np
import math

def overlap(y1,y2,dx):
    v=0.0
    for n in range(len(y1)):
        v+=dx*y1[n]*y2[n]
    return v

def fcf(y1,y2):
    v=0.0;
    norm1=math.sqrt(overlap(y1,y1,6.0/100.0))
    norm2=math.sqrt(overlap(y2,y2,6.0/100.0))
    
    return abs((1.0/norm1)*(1.0/norm2)*overlap(y1,y2,6.0/100.0))

def getallfcf(f1, f2,n1,n2):
    file1=open(f1)
    file2=open(f2)

    lines1=file1.readlines()
    lines2=file2.readlines()
    ys=[[],[]]
    for n in range(n1):
        y=[]
        for i in range(100):
            y.append(float(lines1[i+2].split()[n]))
        ys[0].append(y)
    for n in range(n2):
        y=[]
        for i in range(100):
            y.append(float(lines2[i+2].split()[n]))
        ys[1].append(y)

    
    fcfs=[]
    for i in range(n1):
        fcf_temp=[]
        for j in range(n2):
            fcf_temp.append(fcf(ys[0][i],ys[1][j]))
        fcfs.append(fcf_temp)

    return fcfs

def simulate(fcfs):
    #graph("../Output/Morse1",4,-.5,10,False)
    #graph("../Output/Morse2",4,-.5,10,False)
    f1=open("../Output/Morse1")
    f2=open("../Output/Morse2")
    lines1=f1.readlines()
    lines2=f2.readlines()

    e1=lines1[1].split()
    e2=lines2[1].split()
    
    ax=plt.subplot(111)
    x=[0,0]
    y=[0,1]
    dist=[.9,.1]
    for n in range(2):
        for i in range(10):
            ediff=float(e2[i])-float(e1[n])
            x[0]=ediff
            x[1]=ediff
            y[1]=fcfs[n][i]*dist[n]
            ax.plot(x,y)
    ax.spines[['top','right']].set_visible(False)
    plt.xlabel("E / Eh")
    plt.ylabel("Intensity")
    plt.xlim(0,10)
    plt.show()




def graph(filename, n,ymin,ymax,shouldshow=True):
    f=open(filename)
    
    lines=f.readlines()

    vs=lines[0].split()
    evals=lines[1].split()

    x=np.linspace(-5,5,100);

    ax=plt.subplot(111)
    V=[]
    for i in range(100):
        V.append(float(vs[i]))
    ax.plot(x,V)

    for n in range(n):
        print(evals[n])

    for s in range(n):
        y=[]
        for i in range(100):
            y.append(float(lines[i+2].split()[s])+float(evals[s]))

        ax.plot(x,y)
        plt.hlines(float(evals[s]),-5,5,linestyles='dashed')
    plt.ylim(ymin,ymax)

    plt.xlabel("x / AU")
    plt.ylabel("E / Eh")

    ax.spines[['top','right']].set_visible(False)
    if(shouldshow):
        plt.show()

def main():
    #graph("../Output/Pib", 6,-.2,2)
    #graph("../Output/Fin", 2,-.2,3)
    #graph("../Output/Pbrec", 7,-.2,3.5)
    #graph("../Output/Harm", 7,-.2,5)
    #graph("../Output/Morse1",4,-.5,10,False)
    #graph("../Output/Morse2",4,-.5,10)

    fcfs=getallfcf("../Output/Morse1","../Output/Morse2",2,10)
    #for n in range(len(fcfs)):
    #    print("FRANK CONDON FACTORS ", n)
    #    for m in range(len(fcfs[n])):
    #        print(fcfs[n][m])
    simulate(fcfs)

main()

