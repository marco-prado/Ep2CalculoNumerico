#imports das bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def rk4(x0, n, I, f):
    t = []
    x = [x0]
    h = (I[1] - I[0]) / n
    for k in range(n + 1):
        t.append(I[0] + h * k)
    for k in range(1, n+1):
        K1 = f(t[k - 1], x[k - 1])
        K2 = f(t[k - 1] + h / 2, x[k - 1] + K1 * h / 2)
        K3 = f(t[k - 1] + h / 2, x[k - 1] + K2 * h / 2)
        K4 = f(t[k - 1] + h, x[k - 1] + K3 * h)
        x.append(x[k - 1] + (K1 + 2 * K2 + 2 * K3 + K4) * h / 6)
    return x

def solExplicita(I, n):
    xexplicito = []
    t = 0
    h = (I[1] - I[0]) / n
    while (t <= I[1]+h):
        emenost = np.exp(-t)
        emenos3t = np.exp(-3 * t)
        sint = np.sin(t)
        cost = np.cos(t)
        sin3t = np.sin(3 * t)
        cos3t = np.cos(3 * t)
        xexplicito.append([(emenost * sint) + (emenos3t * cos3t), (emenost * cost) + (emenos3t * sin3t), (-1 * emenost * sint) + (emenos3t * cos3t), (-1 * emenost * cost) + (emenos3t * sin3t)])
        t += h
    return xexplicito

def errork4(xrk4, xexplicito):
    erro = []
    for c in range(len(xrk4)):
        errotemp = np.abs(xexplicito[c][0]-xrk4[c][0])
        for i in range(len(xexplicito[0])):
            if(np.abs(xexplicito[c][i] - xrk4[c][i]) > errotemp):
                errotemp = np.abs(xexplicito[c][i]-xrk4[c][i])
        erro.append(errotemp)
    return erro


#calculo rk4
x0 = [1,1,1,-1]
I = [0, 2]
n = int(input('Insira o valor de n: '))
matA = [[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]]
f = lambda t, x: np.dot(matA,x)

xrk4 = rk4(x0, n, I, f)
print('Solução calculada:')
print(xrk4[len(xrk4) - 1])
print('----------------------------------------------------------')


#calculo solucao explicita
xexplicito = solExplicita(I, n)
print('Solução explícita:')
print(xexplicito[len(xexplicito)-1])
print('----------------------------------------------------------')

#calculo Rs
erros = []
erros.append(max(errork4(rk4(x0, 20, I, f), solExplicita(I,20))))
erros.append(max(errork4(rk4(x0, 40, I, f), solExplicita(I,40))))
erros.append(max(errork4(rk4(x0, 80, I, f), solExplicita(I,80))))
erros.append(max(errork4(rk4(x0, 160, I, f), solExplicita(I,160))))
erros.append(max(errork4(rk4(x0, 320, I, f), solExplicita(I,320))))
erros.append(max(errork4(rk4(x0, 640, I, f), solExplicita(I,640))))

Rs = []
i = 0
while(i<5):
    Rs.append(erros[i]/erros[i+1])
    i += 1

print('Rs calculados:')
print(Rs)

#calculo do erro
t = []
h = (I[1] - I[0]) / n
for k in range(n + 1):
    t.append(I[0] + h * k)

erro = errork4(xrk4, xexplicito)
plt.plot(erro)
plt.title("Gráfico dos Erros para n = " + str(n))
plt.show()

