# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import matplotlib.pyplot as plt
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

# Atribuicoes iniciais
n = 5000
T = [0, 10.0]
h = (T[1] - T[0])/n
f = lambda x, y: (2/3)*x - (4/3)*x*y
g = lambda x, y: (x-1)*y
x0 = 1.5
y0 = 1.5

x, y = eulerExplicito(x0, y0, h, f, g)

t_eixo = []
for i in range(0, n+1):
    t_eixo.append(T[0]+h*i)

plt.subplot(1, 2, 1)
plt.plot(x, y)
plt.title("Retrato de Fase\n")

plt.subplot(1, 2, 2)
plt.plot(t_eixo, x, color='blue')
plt.plot(t_eixo, y, color='red')
plt.title("População x Tempo\n")

plt.show()
