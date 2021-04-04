# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import matplotlib.pyplot as plt
import numpy as np


def eulerExplicito(x0, y0, h, f, g):
    # x0 e y0: valores iniciais
    # h: tamanho do "passo" de cada iteracao
    # f e g = Funções do Sistema EDO

    x = [x0]
    y = [y0]

    for l in range(0, n):
        x.append(x[l] + f(x[l], y[l]) * h)  # x(l+1) = x(l) + f(x(l), y(l))*h
        y.append(y[l] + g(x[l], y[l]) * h)  # y(l+1) = y(l) + g(x(l), y(l))*h

    return x, y


def rk4(x0, y0, h, f, g):
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


# Atribuicoes iniciais
n = 5000
T = [0, 10.0]
h = (T[1] - T[0]) / n
f = lambda x, y: (2 / 3) * x - (4 / 3) * x * y
g = lambda x, y: (x - 1) * y
x0 = 1.5
y0 = 1.5
metodo = int(input("Selecione o método a ser utilizado: \n 1: Método de Euler \n 2: Método de Runge-Kutta de 4a Ordem\n"))

if(metodo == 1):
    x, y = eulerExplicito(x0, y0, h, f, g)
elif(metodo == 2):
    x, y = rk4(x0, y0, h, f, g)

t_eixo = []
for i in range(0, n + 1):
    t_eixo.append(T[0] + h * i)


grafs, (g1, g2) = plt.subplots(1, 2, constrained_layout=True, sharey=True)
g1.plot(x, y)
g1.set_title("Retrato de Fase")

g2.plot(t_eixo, x, color='blue')
g2.plot(t_eixo, y, color='red')
g2.set_title("População x Tempo")

if metodo == 1:
    grafs.suptitle("Método de Euler", fontsize = 15)
elif metodo == 2:
    grafs.suptitle("Método de Runge-Kutta de 4a Ordem", fontsize = 17)

plt.show()