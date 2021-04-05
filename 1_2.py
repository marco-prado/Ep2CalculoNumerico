# MAP 3122 - EP 2 - Tomografia Computadorizada - Métodos numéricos para resolução de EDOs
# Marco Aurélio C. O. Prado - NUSP 11257605
# Silas Lima e Silva        - NUSP 11262233

from math import fabs
import matplotlib.pyplot as plt


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

# Condicoes do Enunciado
x0 = -8.79
T = [1.1, 3.0]
N = 5000
f = lambda t, x: 2*t + (x-t**2)**2

x_num = eulerImplicito(x0, f, N, T)
x_expl = []
E2 = []

tempos = []
h = (T[1]-T[0])/N
for i in range(0, N+1):
    tempos.append(T[0]+h*i)

for i in range(0, N+1):
    x_expl.append(tempos[i]**2+(1/(1-tempos[i])))  # Calculo da solucao explicita
    E2.append(fabs(x_expl[i]-x_num[i]))  # Calculo do Erro na posicao i

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
