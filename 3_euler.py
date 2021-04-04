# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import matplotlib.pyplot as plt


def eulerExplicito(x0, y0, z0, h, f, g, u):
    # x0, y0 e z0: valores iniciais
    # h: tamanho do "passo" de cada iteracao
    # f, g e u = Funções do Sistema EDO

    x = [x0]
    y = [y0]
    z = [z0]

    for l in range(0, n):
        x.append(x[l] + f(x[l], y[l], z[l])*h)  # x(l+1) = x(l) + fx(x(l), y(l), z(l))*h
        y.append(y[l] + g(x[l], y[l], z[l])*h)  # y(l+1) = y(l) + fy(x(l), y(l), z(l))*h
        z.append(z[l] + u(x[l], y[l], z[l])*h)  # z(l+1) = z(l) + fz(x(l), y(l), z(l))*h

    return [x, y, z]

# Atribuicoes iniciais
x0 = 500  # Coelhos
y0 = 500  # Lebres
z0 = 10   # Raposas

n = 5000
T = [0, 100]
alfas = [0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055]
Tfs = [100, 100, 500, 500, 2000, 2000]
retratos_fase = []

for i in range(0, len(alfas)):
    T[1] = Tfs[i]
    alfa = alfas[i]

    h = (T[1] - T[0])/n

    f = lambda x, y, z: x*(1-0.001*x-0.001*y-0.015*z)
    g = lambda x, y, z: y*(1-0.0015*x-0.001*y-0.001*z)
    u = lambda x, y, z: z*(-1+alfa*x+0.0005*y)

    retratos_fase.append(eulerExplicito(x0, y0, z0, h, f, g, u))
