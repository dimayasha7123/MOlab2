from scipy.optimize import golden
from tabulate import tabulate
from math import sqrt
from sympy import *
from enum import Enum


class Type(Enum):
    FAST = 1
    RIVS = 2
    COORD = 3


def countAntiGrad(X):
    return [eval(str(-1 * i)) for i in
            [diff(f, x).subs([(x, X[0]), (y, X[1])]), diff(f, y).subs([(x, X[0]), (y, X[1])])]]


def method(X, eps, type):
    k = 1
    antiGrad = countAntiGrad(X)
    table = [['N', 'x(n)', 'S(n)'], [k, f'({X[0]};{X[1]})', f'({antiGrad[0]};{antiGrad[1]})']]
    while sqrt(sum([pow(i, 2) for i in antiGrad])) >= eps:
        min = golden(lambda z: lambdify((x, y), f)(X[0] + antiGrad[0] * z, X[1] + antiGrad[1] * z))
        for i in range(len(X)):
            X[i] += antiGrad[i] * min
        if type == Type.FAST:
            antiGrad = countAntiGrad(X)
        elif type == Type.RIVS:
            newAntiGrad = countAntiGrad(X)
            B = pow(sum([pow(i, 2) for i in newAntiGrad]), 2) / pow(sum([pow(i, 2) for i in antiGrad]), 2)
            for i in range(len(antiGrad)):
                antiGrad[i] = newAntiGrad[i] + B * antiGrad[i]
        else:
            print("Пока не сделано")
        k += 1
        table.append([k, f'({X[0]};{X[1]})', f'({antiGrad[0]};{antiGrad[1]})'])

    print(tabulate(table, tablefmt='pipe', stralign='center', headers='firstrow'))
    print(f'Всего итераций: {k}')


x, y = symbols('x,y')
f = 4 * pow(x, 2) + 3 * pow(y, 2) - 4 * x * y + x

eps = pow(10, -3)

print('Метод наискорешйего спуска:')
method([0, 0], eps, Type.FAST)

print('Метод сопряженных градиентов (Флетчера-Ривса:')
method([0, 0], eps, Type.RIVS)

# Ответ на стр 92
# x1 = -3/16 = -0.1875
# x2 = -1/8 = -0.125

#изменить колонки таблицы
#изменить визуал таблицы
#обновление в Ривсе каждые k раз
#прикрутить время работы

