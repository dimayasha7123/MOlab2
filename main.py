from scipy.optimize import golden
from tabulate import tabulate
from math import sqrt
from sympy import *

x, y = symbols('x,y')
f = pow(x, 2) + pow((y - 3), 2)


def getAgrad(x1, x2):
    return [agrad(diff(f, x), x1, x2), agrad(diff(f, y), x1, x2)]


def processing(v, eps, x1, x2):  # v - выбор методы : 1 метод наискорейшего спуска, 2 метод сопреженных гардиентов
    k, b1, b2 = 0, 0, 0
    agradx1, agradx2 = getAgrad(x1, x2)
    table = [['N', 'x(n)', 'S(n)'], [k, f'({x1};{x2})', f'({agradx1};{agradx2})']]
    while sqrt(agradx1 ** 2 + agradx2 ** 2) >= eps:
        k += 1
        min = golden(lambda z: f.subs([(x, x1 + agradx1 * z), (y, x2 + agradx2 * z)]))
        x1 += agradx1 * min
        x2 += agradx2 * min
        if v == 1:
            agradx1, agradx2 = getAgrad(x1, x2)
        else:
            b1, b2 = pow(getAgrad(x1, x2)[0], 2) / pow(agradx1, 2), pow(getAgrad(x1, x2)[1], 2) / pow(agradx2, 2)
            agradx1, agradx2 = getAgrad(x1, x2)[0] + b1 * agradx1, getAgrad(x1, x2)[1] + b2 * agradx2
        table.append([k, f'({x1};{x2})', f'({agradx1};{agradx2})'])
        print("a")

    print(sqrt(pow(agradx1, 2) + pow(agradx2, 2)))
    print(tabulate(table, tablefmt='pipe', stralign='center', headers='firstrow'))


def agrad(x, x1, x2):
    return -1 * eval(str(x).replace(" ", "").replace("x", str(x1)).replace("y", str(x2)))


processing(1, 0.0000001, 20, -34), processing(2, 0.0000001, 20, -34)
