# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import matplotlib.pyplot as plt
import numpy as np


def eulerExplicito(x0, y0, n, h, f, g):
    # x0 e y0: valores iniciais
    # h: tamanho do "passo" de cada iteracao
    # f e g = Funções do Sistema EDO

    x = [x0]
    y = [y0]

    for l in range(0, n):
        x.append(x[l] + f(x[l], y[l]) * h)  # x(l+1) = x(l) + f(x(l), y(l))*h
        y.append(y[l] + g(x[l], y[l]) * h)  # y(l+1) = y(l) + g(x(l), y(l))*h

    return x, y


def rk4(x0, y0, n, h, f, g):
    x = [x0]  # array para os valores de x
    y = [y0]  # array para os valores de y
    for k in range(1, n + 1):  # calculo das solucoes
        K1x = f(x[k - 1], y[k - 1])  # calcula K1x
        K1y = g(x[k - 1], y[k - 1])  # calcula K1y
        K2x = f(x[k - 1] + K1x * h / 2, y[k - 1] + K1y * h / 2)  # calcula K2x
        K2y = g(x[k - 1] + K1x * h / 2, y[k - 1] + K1y * h / 2)  # calcula K2y
        K3x = f(x[k - 1] + K2x * h / 2, y[k - 1] + K2y * h / 2)  # calcula K3x
        K3y = g(x[k - 1] + K2x * h / 2, y[k - 1] + K2y * h / 2)  # calcula K3y
        K4x = f(x[k - 1] + K3x * h, y[k - 1] + K3y * h)  # calcula K4x
        K4y = g(x[k - 1] + K3x * h, y[k - 1] + K3y * h)  # calcula K4y
        x.append(x[k - 1] + (K1x + 2 * K2x + 2 * K3x + K4x) * h / 6)  # calcula x e o adiciona ao array
        y.append(y[k - 1] + (K1y + 2 * K2y + 2 * K3y + K4y) * h / 6)  # calcula y e o adiciona ao array
    return x, y  # retorna array calculado

def eulerImplicito(x0, y0, f, g, n, T):
    # x0: valor inicial de x
    # y0: valor inicial de y
    # f e g: Funcoes da EDO para x' e y'
    # n : numero de passos
    # T = [T0,Tf]: Instantes inicial e final

    x = [x0]
    y = [y0]
    h = (T[1]-T[0])/n

    l = 0
    t = T[0]
    for i in range(n):
        xl1_aprox = x[l] + f(x[l], y[l])*h  # Calcula a aproximacao inicial para X_l+1
        yl1_aprox = y[l] + g(x[l], y[l])*h  # Calcula a aproximacao inicial para Y_l+1

        Jacobiano = [[1 - h*(2/3 - (4/3)*y[l]), -h*(-(4/3)*x[l])], [-h*y[l], 1 - h*(x[l] - 1)]]
        invJ = np.linalg.inv(Jacobiano)
        for i in range(0, 7):  # 7 passos do metodo de newton
            xl1_new = xl1_aprox - invJ[0][0]*(xl1_aprox-h*f(xl1_aprox, yl1_aprox) - x[l]) - invJ[0][1]*(yl1_aprox-h*g(xl1_aprox, yl1_aprox) - y[l])
            yl1_aprox = yl1_aprox - invJ[1][0]*(xl1_aprox-h*f(xl1_aprox, yl1_aprox) - x[l]) - invJ[1][1]*(yl1_aprox-h*g(xl1_aprox, yl1_aprox) - y[l])
            xl1_aprox = xl1_new
        x.append(x[l] + f(xl1_aprox, yl1_aprox)*h)  # Implementa o metodo de Euler implicito: X_l+1 = X_l + f(t_l+1, x_l+1, y_l+1)*h
        y.append(y[l] + g(xl1_aprox, yl1_aprox)*h)  # Implementa o metodo de Euler implicito: Y_l+1 = Y_l + f(t_l+1, x_l+1, y_l+1)*h
        l += 1
    return x, y


# Atribuicoes iniciais
n = 5000
T = [0, 10.0]
h = (T[1] - T[0]) / n
f = lambda x, y: (2 / 3) * x - (4 / 3) * x * y
g = lambda x, y: (x - 1) * y
x0 = 1.5
y0 = 1.5

exercicio = int(input("Selecione o exercicio: \n 1: Resolução Utilizando o Método de Euler Explícito \n 2: Resolução Utilizando o Método de Euler Implícito \n 3: Comparação Euler Explícito x Implícito \n 4: Resolução Utilizando o Método de Runge-Kutta de 4a Ordem\n"))

if(exercicio == 1 or exercicio == 2 or exercicio == 4):
    if(exercicio == 1):
        x, y = eulerExplicito(x0, y0, n, h, f, g)
    elif (exercicio == 2):
        n = 750
        x, y = eulerImplicito(x0, y0, f, g, n, T)
    elif(exercicio == 4):
        x, y = rk4(x0, y0, n, h, f, g)

    t_eixo = []
    for i in range(0, n + 1):
        t_eixo.append(T[0] + h * i)

    grafs, (g1, g2) = plt.subplots(1, 2, constrained_layout=True, sharey=True)
    g1.plot(x, y)
    g1.set_title("Retrato de Fase")

    g2.plot(t_eixo, x, color='blue')
    g2.plot(t_eixo, y, color='red')
    g2.set_title("População x Tempo")

    if exercicio == 1:
        grafs.suptitle("Método de Euler Explícito", fontsize = 15)
    elif exercicio == 2:
        grafs.suptitle("Método de Euler Implícito", fontsize = 17)
    elif exercicio == 4:
        grafs.suptitle("Método de Runge-Kutta de 4a Ordem", fontsize = 17)

    plt.show()