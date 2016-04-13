from __future__ import print_function
import math
import matplotlib.pyplot as plt

test = 1  # set to 0 to switch off

def antoine (length):
    p = list()
    for l in length:
        p.append(1.0/(400**l))
    return(p)

def same_amino_acid_sequence (length):
    p = list()
    for l in length:
        p.append((1.0/20)**l)
    return(p)

def two_same_birthday (x, n):
    '''
    Description: calculates probability of two people having the same birthday
    In: x is list of values (for example nr of days in a year), n = nr of people
    In x is list of options for CDR3 sequence, n = nr of vgenes
    Out: list of probabilities
    NOTE: dit geeft een OverflowError
    '''
    p = list()
    for m in x:
        m = float(m)
        p.append(1 - (math.factorial(m) / (m**n * math.factorial(m-n))))
    return(p)

def taylor_approximation(x, n):
    '''
    Description: calculates probability of two people having the same birthday
    In: x is list of values (for example nr of days in a year), n = nr of people
    In x is list of options for CDR3 sequence, n = nr of vgenes
    Out: list of probabilities
    '''
    p = list()
    for m in x:
        m = float(m)
        p.append(1-math.exp((-1*n**2)/(2*m)))
    return(p)

def same_as_you (x, n):
    '''
    Description: calculates probability of someone having same birthday, given a certain date
    In: x is list of values (for example nr of days in a year), n = nr of people
    In x is list of options for CDR3 sequence, n = nr of vgenes
    Out: list of probabilities
    '''
    p = list()
    for m in x:
        m = float(m)
        p.append(1.0 - ((m-1)/m)**n)
    return(p)

def autolabel(x,y):
    # attach some text labels
    for i in range(len(x)):
        ax.text(x[i], y[i],'%d' % y[i],ha='center', va='bottom')

if test == 1:

    # Voorbeeld wikipedia (Birthday Problem)
    # m = 365.0
    # kids = range(1,367)
    # p = list()
    # for n in kids:
    #     n = float(n)
    #     p.append(1.0 - ((m-1)/m)**n)
    # x = kids
    # y = p

    x = range(1,11)
    cdr3s = [20**a for a in x]
    nr_vgenes = 65                  # TCRb: 52, BCRh: 65
    y = same_as_you(cdr3s, nr_vgenes)
    z = taylor_approximation(cdr3s, nr_vgenes)
    w = same_amino_acid_sequence(x)
    v = antoine(x)
    print(x)
    print(y)
    print(z)
    print(w)
    print(v)

    fig, ax = plt.subplots()
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    p1 = ax.plot(x, y, color="blue", linestyle="solid", marker="o")
    p2 = ax.plot(x, z, color="red", linestyle="dashed", marker="o")
    p3 = ax.plot(x, w, color="green", linestyle="dashdot", marker="o")
    p4 = ax.plot(x, v, color="purple", linestyle="dotted", marker="o")
    ax.set_xlabel("Number of inserted amino acids in CDR3")
    ax.set_ylabel("Probability")
    ax.set_title("Probability of two V genes having the same CDR3")
    plt.xticks(x)
    plt.legend(["Same V gene given a CDR3","taylor_approximation","(1/20)^L","1.0/(400^l)"])
    plt.show()
