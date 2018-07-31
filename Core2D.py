import numpy as np
import matplotlib.pyplot as plt
import gc

d=1.0
L=41.0
r=0.0
vac=0.0

##N=int((L+2*r)/d)
n_vac=int(vac/d)
n_r=int(r/d)
n_c=int(L/d)
n=2*n_vac+2*n_r+n_c
print n_vac
print n_r
print n_c
print n
N=(n+2)**2

Dvac=0.16
Svac=1e8

##Dr=0.16
##Sar=0.02
Dr=0.16
Sar=0.02
Dp=0.16
Sap=20.0
Dc=9.21
Sac=0.1532
nuSf_val=0.157

######################################################      
##'Utworzenie macierzy nuSf (odpowiadajacej ukladowi X Y)
##'i przerobienie na wektor
nuSf_mac=np.zeros(((n+2),(n+2)))

for i in range (1,n+1):
    if i>(n_vac+n_r) and i<=(n_vac+n_r+n_c):
        nuSf_mac[i,n_vac+n_r+1:(n_vac+n_r+n_c)+1]=nuSf_val
##    nuSf_mac[i,i]=1
##    nuSf_mac[i,n+1-i]=1
nuSf=np.zeros(n**2)
for j in range (1,(n+1)):
    for i in range (1,(n+1)):
        nuSf[i-1+n*(j-1)]=nuSf_mac[j,i]
######################################################      
##'Utworzenie macierzy D (odpowiadajacej ukladowi X Y)
##'i przerobienie na wektor
D_mac=np.zeros(((n+2),(n+2)))
D_mac[0,:]=Dvac
D_mac[(n+2)-1,:]=Dvac
for i in range (1,n+1):       
    D_mac[i,:]=Dvac   
    if i>n_vac and i<=(n_vac+2*n_r+n_c):
        D_mac[i,n_vac+1:n+1-n_vac]=Dr
    if i>(n_r+n_vac) and i<=(n_r+n_vac+n_c):
        D_mac[i,n_r+n_vac+1:(n_r+n_vac+n_c)+1]=Dc
##    D_mac[i,i]=1
##    D_mac[i,n+1-i]=1
######################################################      
##'Utworzenie macierzy Sa (odpowiadajacej ukladowi X Y)
##'i przerobienie na wektor
Sa_mac=np.zeros(((n+2),(n+2)))
Sa_mac[0,:]=Svac
Sa_mac[(n+2)-1,:]=Svac
for i in range (1,n+1):       
    Sa_mac[i,:]=Svac      
    if i>n_vac and i<=(n_vac+2*n_r+n_c):
        Sa_mac[i,n_vac+1:n+1-n_vac]=Sar
    if i>(n_vac+n_r) and  i<=(n_vac+n_r+n_c):
        Sa_mac[i,n_vac+n_r+1:(n_vac+n_r+n_c)+1]=Sac
##    Sa_mac[i,i]=0.5
##    Sa_mac[i,n+1-i]=0.5
######################################################
##D_mac2=[D,D]
##Sa_mac2=[Sa,Sa,Sa,Sa,Sa]
##nuSf_mac2=[nuSf,nuSf]
#####wyrysowanie macierzy
##plt.matshow(D_mac)
##plt.matshow(D_mac2)
##plt.matshow(Sa_mac)
##plt.matshow(Sa_mac2)
##plt.matshow(nuSf_mac)
##plt.matshow(nuSf_mac2)
##plt.show()
######################


for y in range (0,n,4):
    ########PRET KONTROLNY##########
    ###glebokosc preta:
##    y=
    ###polozenie preta
    x=int(round(n/2.0))
    ###wypelnienie macierzy
    for i in range(1+n_vac,y+1+n_vac):
        D_mac[i,x]=Dp
        Sa_mac[i,x]=Sap
    ##plt.matshow(Sa_mac)
    ##plt.show()


    ran=range(1+(n+2),(n**2+3*n)+1)
    N_B=n**2+2*n-2
    ran_del=range(0,n+3)
    ran_A=range(0,n**2)


    F=np.ones(n**2, dtype=np.float64)
    S=np.zeros(n**2, dtype=np.float64)
    A=np.zeros((n**2,n**2), dtype=np.float64)

    for k in ran_A:
        i=(k%n)+1
        j=int(k/n)+1
    ##    print k
    ##    print i
    ##    print j
        if k-n>=0:
            A[k,k-(n)]=-1/(2*d**2)*(D_mac[j,i]+D_mac[j-1,i])
        if k-1>=0:
            A[k,k-1] = -1/(2*d**2)*(D_mac[j,i]+D_mac[j,i-1])
        A[k,k] = 1/(2*d**2)*(D_mac[j+1,i]+D_mac[j,i+1]+4*D_mac[j,i]
                             +D_mac[j,i-1]+D_mac[j-1,i])+Sa_mac[j,i]
        if k+1<n**2:
            A[k,k+1] = -1/(2*d**2)*(D_mac[j,i]+D_mac[j,i+1])
        if k+n<n**2:
            A[k,k+(n)]=-1/(2*d**2)*(D_mac[j,i]+D_mac[j+1,i])
##    print A
    #wyrysowanie macierzy)
    ##plt.matshow(A)
    ##plt.show()
    ######################
##    del(nuSf_mac)
##    del(D_mac)
##    del(Sa_mac)

    gc.collect

    lmd = 1.0
    S=1.0/lmd*nuSf*F
    Fnew=np.copy(F)
    eps=1e-6

    while True:
    ##    print Fnew[0:n-1]
    ##    for i in range(0,n):
    ##        print Fnew[i*n]
        Fnew=np.linalg.solve(A,S)
        lmd_new=sum(nuSf*Fnew*d)/sum(nuSf*F)*lmd
        S=nuSf*Fnew/lmd_new
        F=np.copy(Fnew)
         
        if abs((lmd_new-lmd)/lmd_new)<eps:
            break
        lmd=lmd_new
        print lmd
    print 'Policzono'
    print int(y*100/n+1)
    F_mac=np.zeros((n,n))
    for k in range(0,n**2):
        F_mac[int(k/(n)),k%(n)]=F[k]

#    plt.plot(F_mac[int((n)/2),0:n])
    plt.plot(F_mac[0:n,int(round((n)/2.0))],'--')
    plt.ylabel('Y Flux profile')
    plt.xlabel('X [cm]')


#    flux=np.sum(F_mac)
#    plt.plot(int(y*100/n+1),flux, 'bo')
#    plt.ylabel('Total flux')
#    plt.xlabel('Rod position')

##    plt.xlim( 0, 100 )


plt.show()

