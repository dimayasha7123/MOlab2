from scipy.optimize import golden
from tabulate import tabulate
from math import sqrt
from sympy import *
import time


def countAntiGrad(X):
    return [eval(str(-1 * i)) for i in
            [diff(f, x).subs([(x, X[0]), (y, X[1])]), diff(f, y).subs([(x, X[0]), (y, X[1])])]]


def method(X0, eps, type, n=3):
    X = X0.copy()
    startTime = time.time()
    k, B = 0, None
    antiGrad = countAntiGrad(X)
    table = [['k', 'X1(n)', 'X2(n)', 'S(n)'], [k, f'{X[0]}', f'{X[1]}', f'({antiGrad[0]}; {antiGrad[1]})']]
    if type == 'rivs':
        table = [['k', 'X1(n)', 'X2(n)', 'S(n)', 'B'],
                 [k, f'{X[0]}', f'{X[1]}', f'({antiGrad[0]}; {antiGrad[1]})', '-']]
    while sqrt(sum([pow(i, 2) for i in antiGrad])) >= eps:
        k += 1
        min = golden(lambda z: lambdify((x, y), f)(X[0] + antiGrad[0] * z, X[1] + antiGrad[1] * z))
        for i in range(len(X)):
            X[i] += antiGrad[i] * min
        if type == 'fast':
            antiGrad = countAntiGrad(X)
        elif type == 'rivs':
            newAntiGrad = countAntiGrad(X)
            B = pow(sum([pow(i, 2) for i in newAntiGrad]), 2) / pow(sum([pow(i, 2) for i in antiGrad]), 2)
            if (k + 1) % n == 0:
                B = 0
            for i in range(len(antiGrad)):
                antiGrad[i] = newAntiGrad[i] + B * antiGrad[i]
        else:
            print("Пока не сделано")
        if type == 'fast':
            table.append([k, f'{X[0]}', f'{X[1]}', f'({antiGrad[0]};{antiGrad[1]})'])
        elif type == 'rivs':
            table.append([k, f'{X[0]}', f'{X[1]}', f'({antiGrad[0]};{antiGrad[1]})', B])
        else:
            table.append(['Not exist'])

    print(tabulate(table, tablefmt='fancy_grid', numalign='center', stralign='center', headers='firstrow', floatfmt=(".1d", ".10f", ".10f")))
    print(f'Всего итераций: {k + 1}')
    print(f'Время выполнения: {(time.time() - startTime)}')


if __name__ == '__main__':
    x, y = symbols('x,y')
    f = 4 * pow(x, 2) + 3 * pow(y, 2) - 4 * x * y + x

    Xinput = [0, 0]
    eps = pow(10, -7)

    print('Метод наискорешйего спуска:')
    method(Xinput, eps, 'fast')

    print('Метод сопряженных градиентов (Флетчера-Ривса):')
    method(Xinput, eps, 'rivs', 4)

    # Ответ на стр 92
    # x1 = -3/16 = -0.1875
    # x2 = -1/8 = -0.125


    # копировать X на входе в функцию
