# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

import matplotlib.pyplot as plt
import numpy as np

def eulerExplicito(x0, y0, z0, n, h, f, g, u):
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

def rk4(x0, y0, z0, n, h, f, g, u):
    x = [x0]  # array para os valores de x
    y = [y0]  # array para os valores de y
    z = [z0]  # array para os valores de z
    for k in range(1, n + 1):  # calculo das solucoes
        K1x = f(x[k - 1], y[k - 1], z[k - 1])  # calcula K1x
        K1y = g(x[k - 1], y[k - 1], z[k - 1])  # calcula K1y
        K1z = u(x[k - 1], y[k - 1], z[k - 1])  # calcula K1z
        K2x = f(x[k - 1] + K1x * h / 2, y[k - 1] + K1y * h / 2, z[k - 1] + K1z * h / 2)  # calcula K2x
        K2y = g(x[k - 1] + K1x * h / 2, y[k - 1] + K1y * h / 2, z[k - 1] + K1z * h / 2)  # calcula K2y
        K2z = u(x[k - 1] + K1x * h / 2, y[k - 1] + K1y * h / 2, z[k - 1] + K1z * h / 2)  # calcula K2z
        K3x = f(x[k - 1] + K2x * h / 2, y[k - 1] + K2y * h / 2, z[k - 1] + K2z * h / 2)  # calcula K3x
        K3y = g(x[k - 1] + K2x * h / 2, y[k - 1] + K2y * h / 2, z[k - 1] + K2z * h / 2)  # calcula K3y
        K3z = u(x[k - 1] + K2x * h / 2, y[k - 1] + K2y * h / 2, z[k - 1] + K2z * h / 2)  # calcula K3z
        K4x = f(x[k - 1] + K3x * h, y[k - 1] + K3y * h, z[k - 1] + K3z * h)  # calcula K4x
        K4y = g(x[k - 1] + K3x * h, y[k - 1] + K3y * h, z[k - 1] + K3z * h)  # calcula K4y
        K4z = u(x[k - 1] + K3x * h, y[k - 1] + K3y * h, z[k - 1] + K3z * h)  # calcula K4y

        x.append(x[k - 1] + (K1x + 2 * K2x + 2 * K3x + K4x) * h / 6)  # calcula x e o adiciona ao array
        y.append(y[k - 1] + (K1y + 2 * K2y + 2 * K3y + K4y) * h / 6)  # calcula y e o adiciona ao array
        z.append(z[k - 1] + (K1z + 2 * K2z + 2 * K3z + K4z) * h / 6)  # calcula z e o adiciona ao array

    return [x, y, z]  # retorna array calculado



# Atribuicoes iniciais
x0 = 500  # Coelhos
y0 = 500  # Lebres
z0 = 10   # Raposas

n = 8000
T = [0, 100]
alfas = [0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055]
Tfs = [100, 100, 500, 500, 2000, 2000]
retratos_fase = []

metodo = int(input("Selecione o método a ser utilizado: \n 1: Método de Euler \n 2: Método de Runge-Kutta de 4a Ordem\n"))
alfa_selecionado = int(input("Selecione o alfa a ser utilizado:\n 1: 0.001\n 2: 0.002\n 3: 0.0033\n 4: 0.0036\n 5: 0.005\n 6: 0.0055\n")) - 1

if(input("Você deseja alterar o valor inicial padrão de coelhos (500), lebres (500) e raposas (10)?(s/n) ") == 's'):
    x0 = int(input("Insira o novo valor de coelhos: "))
    y0 = int(input("Insira o novo valor de lebres: "))
    z0 = int(input("Insira o novo valor de raposas: "))

if(input("Você deseja alterar o valor padrão de Tf?(s/n) ") == 's'):
    Tfs[alfa_selecionado] = int(input("Insira o novo valor de Tf: "))

for i in range(0, len(alfas)):
    T[1] = Tfs[i]
    alfa = alfas[i]

    h = (T[1] - T[0])/n

    f = lambda x, y, z: x*(1-0.001*x-0.001*y-0.015*z)
    g = lambda x, y, z: y*(1-0.0015*x-0.001*y-0.001*z)
    u = lambda x, y, z: z*(-1+alfa*x+0.0005*y)

    if metodo == 1:
        retratos_fase.append(eulerExplicito(x0, y0, z0, n, h, f, g, u))
    elif metodo == 2:
        retratos_fase.append(rk4(x0, y0, z0, n, h, f, g, u))

print("Numero de coelhos: " + str(np.floor(retratos_fase[alfa_selecionado][0][-1])))
print("Numero de lebres: " + str(np.floor(retratos_fase[alfa_selecionado][1][-1])))
print("Numero de raposas: " + str(np.floor(retratos_fase[alfa_selecionado][2][-1])))

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot3D(retratos_fase[alfa_selecionado][0], retratos_fase[alfa_selecionado][1], retratos_fase[alfa_selecionado][2], 'gray')
plt.title("Retrato com Alfa = " + str(alfas[alfa_selecionado]))
plt.show()