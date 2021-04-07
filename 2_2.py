# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import numpy as np
import matplotlib.pyplot as plt

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
    while t < (T[1]-h):
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
        t += h
    return x, y

# Atribuicoes iniciais
n = 750
T = [0, 10.0]
h = (T[1] - T[0])/n
f = lambda x, y: (2/3)*x - (4/3)*x*y
g = lambda x, y: (x-1)*y
x0 = 1.5
y0 = 1.5

x, y = eulerImplicito(x0, y0, f, g, n, T)

t_eixo = []
for i in range(0, n):
    t_eixo.append(T[0]+h*i)

plt.subplot(1, 2, 1)
plt.plot(x, y)
plt.title("Retrato de Fase\n")

plt.subplot(1, 2, 2)
plt.plot(t_eixo, x, color='blue')
plt.plot(t_eixo, y, color='red')
plt.title("População x Tempo\n")

plt.show()