# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

#imports das bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt
from math import fabs


#funcao que calcula a solucao de uma EDO utilizando o metodo de Runge-Kutta de ordem 4
def rk4(x0, n, intervalo, f):
    x = [x0]  # array para os valores de x
    h = (intervalo[1] - intervalo[0]) / n #passo calculado atraves do intervalo dado e de n
    for k in range(1, n+1): #calculo das solucoes
        K1 = f(x[k - 1]) #calcula K1
        K2 = f(x[k - 1] + K1 * h / 2) #calcula K2
        K3 = f(x[k - 1] + K2 * h / 2) #calcula K3
        K4 = f(x[k - 1] + K3 * h) #calcula K4
        x.append(x[k - 1] + (K1 + 2 * K2 + 2 * K3 + K4) * h / 6) #calcula x e o adiciona ao array
    return x #retorna array calculado


#funcao que calcula os valores de x através da solucao explicita
def solExplicita(intervalo, n):
    xexplicito = [] #array para os valores de x
    t = 0
    h = (intervalo[1] - intervalo[0]) / n #calcula o passo atraves do intervalo dado e de n
    while (t <= intervalo[1]+h):
        #calculo dos valores utilizados na solucao explicita
        emenost = np.exp(-t)
        emenos3t = np.exp(-3 * t)
        sint = np.sin(t)
        cost = np.cos(t)
        sin3t = np.sin(3 * t)
        cos3t = np.cos(3 * t)
        #calculo do valor de x e insercao dele no array das solucoes
        xexplicito.append([(emenost * sint) + (emenos3t * cos3t), (emenost * cost) + (emenos3t * sin3t), (-1 * emenost * sint) + (emenos3t * cos3t), (-1 * emenost * cost) + (emenos3t * sin3t)])
        t += h
    return xexplicito #retorna array calculado

#funcao que calcula o erro entre a solucao calculada atraves do RK4 e a explicita
def errork4(xrk4, xexplicito):
    erro = []
    for c in range(len(xrk4)): #for para percorrer linhas
        errotemp = np.abs(xexplicito[c][0]-xrk4[c][0]) #variavel que armazena o maior erro de cada linha
        for i in range(len(xexplicito[0])): #for para percorrer colunas
            if(np.abs(xexplicito[c][i] - xrk4[c][i]) > errotemp): #caso o valor atual seja maior que errotemp, ele se torna o novo errotemp
                errotemp = np.abs(xexplicito[c][i]-xrk4[c][i])
        erro.append(errotemp) #insere o maior erro no array de erros
    return erro #retorna array calculado

def eulerImplicito(x0, f, n, T):
    # x0: valor inicial de x
    # f: Funcao da EDO
    # n : numero de passos
    # T = [T0,Tf]: Instantes inicial e final

    x = [x0]
    h = (T[1]-T[0])/n
    l = 0
    t = T[0]
    while t < (T[1]-h):
        xl1_aprox = x[l] + f(t, x[l])*h  # Calcula a aproximacao inicial para X_l+1
        for i in range(0, 7):  # 7 passos do metodo de newton
            xl1_aprox = xl1_aprox - ((x[l]-xl1_aprox + f(t+h, xl1_aprox)*h) / (2*h*(xl1_aprox - t**2)-1))
        x.append(x[l] + f(t+h, xl1_aprox)*h)  # Implementa o metodo de Euler implicito: X_l+1 = X_l + f(t_l+1, x_l+1)*h
        l += 1
        t += h
    return x


exercicio = int(input("Selecione o exercício: \n1: Teste Runge-Kutta 4 \n2: Teste Euler Implícito\n"))

if exercicio == 1:
    #calculo rk4
    x0 = [1,1,1,-1]
    intervalo = [0, 2]
    n = int(input('Insira o valor de n: '))
    matA = [[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]]
    f = lambda x: np.dot(matA,x) #f(t, x(t))
    xrk4 = rk4(x0, n, intervalo, f)
    print('Solução Calculada:')
    print(xrk4[len(xrk4) - 1])
    print('----------------------------------------------------------')


    #calculo solucao explicita
    xexplicito = solExplicita(intervalo, n)
    print('Solução Explícita:')
    print(xexplicito[len(xexplicito)-1])
    print('----------------------------------------------------------')


    #calculo Rs
    erros = []
    #calculo do maior erro para cada valor de n
    erros.append(max(errork4(rk4(x0, 20, intervalo, f), solExplicita(intervalo,20))))
    erros.append(max(errork4(rk4(x0, 40, intervalo, f), solExplicita(intervalo,40))))
    erros.append(max(errork4(rk4(x0, 80, intervalo, f), solExplicita(intervalo,80))))
    erros.append(max(errork4(rk4(x0, 160, intervalo, f), solExplicita(intervalo,160))))
    erros.append(max(errork4(rk4(x0, 320, intervalo, f), solExplicita(intervalo,320))))
    erros.append(max(errork4(rk4(x0, 640, intervalo, f), solExplicita(intervalo,640))))
    Rs = []
    i = 0
    while(i<5): #preenchimento do array de Rs
        Rs.append(erros[i]/erros[i+1])
        i += 1
    print('Rs Calculados:')
    print(Rs)


    #calculo e plot do erro
    t = []
    h = (intervalo[1] - intervalo[0]) / n
    for k in range(n + 1):
        t.append(intervalo[0] + h * k)
    erro = errork4(xrk4, xexplicito)
    plt.plot(erro)
    plt.title("Gráfico de erros para n = " + str(n))
    plt.show()

elif exercicio == 2:
    # Condicoes do Enunciado
    x0 = -8.79
    T = [1.1, 3.0]
    N = 5000
    f = lambda t, x: 2 * t + (x - t ** 2) ** 2

    x_num = eulerImplicito(x0, f, N, T)
    x_expl = []
    E2 = []

    tempos = []
    h = (T[1] - T[0]) / N
    for i in range(0, N + 1):
        tempos.append(T[0] + h * i)

    for i in range(0, N + 1):
        x_expl.append(tempos[i] ** 2 + (1 / (1 - tempos[i])))  # Calculo da solucao explicita
        E2.append(fabs(x_expl[i] - x_num[i]))  # Calculo do Erro na posicao i

    plt.subplot(1, 3, 1)
    plt.plot(tempos, x_expl)
    plt.title("Solucao Explicita")

    plt.subplot(1, 3, 2)
    plt.plot(tempos, x_num)
    plt.title("Solucao Numerica")

    plt.subplot(1, 3, 3)
    plt.plot(tempos, E2)
    plt.title("Erros")

    plt.show()
