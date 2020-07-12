# OptimizeCost (M/M/C/K)
# author: Dinh Thi Nhan

import numpy as np
import json

def factorial(x):
    if  x==1 or x==0:
        return 1
    else:
        return x*factorial(x-1)

def Po(K,C,ro):
    Po = 1
    for n in range(1, C):
        Po = Po + ro**n/factorial(n)
    for n in range(C, K+1):
        Po= Po + (C**(n-C)*factorial(C))**(-1)*ro**n
    return  Po**(-1)

def MeanWaitingTime (K,C,ro,lamda):                     # Mean waiting time in queue
    lameff = lamda * (1 - 1/(factorial(C)*(C**(K-C))) * ro**K * Po(K, C, ro))   #lamda effective
    Ro = ro/C
    # Mean numbers of customers in queue
    Lq = Po(K,C,ro)* (C*Ro)**C*Ro/(factorial(C)*(1-Ro)**2) * (1 - Ro**(K-C+1)-(1-Ro)*(K-C+1)*Ro**(K-C))

   # Mean numbers of customers in system
        #L = Lq +C
        #for n in range(C):
            #L = L - (C-n)* ro**n /factorial(n) * Po(K,C,ro)

    return Lq / lameff

def optimize(veclamda,luy,tw,c,Ck, Cc,K_max, C_max):
    lamda = sum(veclamda)
    ro = lamda/luy
    flag = 0
    for C in range(1,C_max):                                # C is numbers of server in system
        for K in range(C+1, K_max):                           # K is maximum numbers of customers in system
            if MeanWaitingTime(K,C,ro,lamda) <= tw :
                Pr =  1 / (factorial(C) * (C ** (K-C))) * ro ** K * Po(K, C, ro)  #Probability a customer is denied
                f_total = Pr* veclamda.dot(np.transpose(c)) + (K-C)*Ck + C*Cc
                flag = 1
            break
        break

    with open('Output.txt', 'w') as f:
        if flag :
            minf = f_total
            K_min, C_min = K,C
            for n in range(C, C_max):
                for m in range(n+1, K_max):
                    if MeanWaitingTime(m, n, ro, lamda) <= tw:
                        Pr = 1 / (factorial(n) * (n ** (m-n))) * ro ** m * Po(m, n, ro)
                        f_total = Pr * veclamda.dot(np.transpose(c)) + (m-n) * Ck + n * Cc
                        if f_total < minf :
                            minf = f_total
                            K_min, C_min = m,n

            text = ("f_total = " + str(minf) + "\nNumbers of servers in system  = " + str(C_min)
                    + '\nMaximum numbers of customers possible in the system = ' + str(K_min)+ '\n')
            f.write(text)
        else:
            f.write('No found K,C within the MeanWaitingTime satisfaction range!\n')


def main():
    with open('Input.txt', 'r') as f :
        json_object = json.load(f)
        Ck            = json_object["Cost per queue"]
        Cc            = json_object["Cost per server"]
        veclamda      = np.array(json_object["veclamda"])   # rate of customers arriving for each customer group
        c             = np.array(json_object["c"])          # cost when a customer in each customers group is denied
        luy           = json_object["Rate of going state n to state n-1"]
        tw            = json_object["Maximum of MeanWaitingTime"]
        K_max         = json_object["K_max"]
        C_max         = json_object["C_max"]
                 # paremeter

    optimize(veclamda, luy, tw, c, Ck, Cc,K_max, C_max)

main()