from scipy.optimize import golden
from tabulate import tabulate
from math import sqrt
from sympy import *
from enum import Enum

x, y = symbols('x,y')
f = pow(x, 2) + pow(y - 3, 2)


def countAntiAgrad(x1, x2):
    print([eval(str(-1 * diff(f, x).subs([(x, x1), (y, x2)]))), eval(str(-1 * diff(f, y).subs([(x, x1), (y, x2)])))])
    return [eval(str(-1 * diff(f, x).subs([(x, x1), (y, x2)]))), eval(str(-1 * diff(f, y).subs([(x, x1), (y, x2)])))]


def processing(x1, x2, eps, type):
    k = 1
    antiGradX1, antiGradX2 = countAntiAgrad(x1, x2)
    table = [['N', 'x(n)', 'S(n)'], [k, f'({x1};{x2})', f'({antiGradX1};{antiGradX2})']]
    while sqrt(pow(antiGradX1, 2) + pow(antiGradX2, 2)) >= eps:
        k += 1
        min = golden(lambda z: lambdify((x, y), f)(x1 + antiGradX1 * z, x2 + antiGradX2 * z))
        x1 += antiGradX1 * min
        x2 += antiGradX2 * min
        if type == Type.FAST:
            antiGradX1, antiGradX2 = countAntiAgrad(x1, x2)
        elif type == Type.RIVS:
            newAntiGradX1, newAntiGradX2 = countAntiAgrad(x1, x2)
            b1, b2 = pow(newAntiGradX1, 2) / pow(antiGradX1, 2), pow(newAntiGradX2, 2) / pow(antiGradX2, 2)
            denAntiGradX1, denAntiGradX2 = newAntiGradX1 + b1 * antiGradX1, newAntiGradX2 + b2 * antiGradX2
            print(f'b1 = {b1}')
            print(f'b2 = {b2}')
            print(f'denAntiGrad = {denAntiGradX1, denAntiGradX2}')


            b = pow(pow(newAntiGradX1, 2) + pow(newAntiGradX2, 2), 2) / pow(pow(antiGradX1, 2) + pow(antiGradX2, 2), 2)
            antiGradX1, antiGradX2 = newAntiGradX1 + b * antiGradX1, newAntiGradX2 + b * antiGradX2
            print(f'b = {b}')
            print(f'antiGrad = {antiGradX1, antiGradX2}')


        table.append([k, f'({x1};{x2})', f'({antiGradX1};{antiGradX2})'])

    print(tabulate(table, tablefmt='pipe', stralign='center', headers='firstrow'))


class Type(Enum):
    FAST = 1
    RIVS = 2
    COORD = 3


a, b = 20, -34
eps = pow(10, -22)
processing(20, -34, eps, Type.FAST)
processing(20, -34, eps, Type.RIVS)