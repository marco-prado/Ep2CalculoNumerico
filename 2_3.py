# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import numpy as np
import matplotlib.pyplot as plt
from math import fabs

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
    for l in range(0, n):
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
        t += h
    return x, y

def eulerExplicito(x0, y0, h, f, g):
# x0 e y0: valores iniciais
# h: tamanho do "passo" de cada iteracao
# f e g = Funções do Sistema EDO

    x = [x0]
    y = [y0]

    for l in range(0, n):
        x.append(x[l] + f(x[l], y[l])*h)  # x(l+1) = x(l) + f(x(l), y(l))*h
        y.append(y[l] + g(x[l], y[l])*h)  # y(l+1) = y(l) + g(x(l), y(l))*h

    return x, y

# Valores fixos
T = [0, 10.0]
f = lambda x, y: (2/3)*x - (4/3)*x*y
g = lambda x, y: (x-1)*y
x0 = 1.5
y0 = 1.5

valoresN = [250, 500, 1000, 2000, 4000]

pos = 1
for n in valoresN:
    xim, yim = eulerImplicito(x0, y0, f, g, n, T)

    h = (T[1] - T[0]) / n
    xex, yex = eulerExplicito(x0, y0, h, f, g)

    Ex = np.subtract(xim, xex)
    Ey = np.subtract(yim, yex)

    t_eixo = []
    for i in range(0, n+1):
        t_eixo.append(T[0] + ((T[1]-T[0])/n)*i)

    plt.subplot(1, 5, pos)
    plt.plot(t_eixo, Ex, color='blue')
    plt.plot(t_eixo, Ey, color='red')
    plt.title("Erro com N = " + str(n))

    pos += 1

plt.show()
