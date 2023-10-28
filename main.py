Задание 1 (ряды Фурье)
import sympy as sp
import matplotlib.pyplot as plt
t = sp.symbols('t')
Num = 21
Gr = 5
f = Gr * (t - Num) ** 2 + t
f = sp.Piecewise((f, (t >= 0) & (t <= 2 * Num)), (0, True))
# Интервал [0, 2*Num]
# Cимвольное разложение функции в ряд Фурье
fourier_series = sp.fourier_series(f, (t, 0, 2 * Num))
# Создание списков для хранения значений t и соответствующих значений f(t),
S0(t), S1(t), S2(t), S10(t)
t_values = []
f_values = []
S_values = {0: [], 1: [], 2: [], 10: []}
# Вычисление значений и заполнение списков
for t_val in range(0, 2 * Num + 1):
 t_values.append(t_val)
 f_val = f.subs(t, t_val)
 f_values.append(f_val)
 for i in [0, 1, 2, 10]:
 S_val = fourier_series.truncate(i).subs(t, t_val)
 S_values[i].append(S_val)
# Построение графиков
for i in [0, 1, 2, 10]:
 plt.figure(figsize=(10, 6))
 plt.plot(t_values, f_values, label='f(t)', linestyle='--', marker='o')
 plt.plot(t_values, S_values[i], label=f'S{i}(t)')
 plt.xlabel('t')
 plt.ylabel('f(t)')
 plt.legend()
 plt.title(f'График функции: f(t) и график частичной суммы: S{i}(t)')
 plt.grid(True)
 plt.show()
# Максимальная погрешность приближения функции f(t) частичной суммой S50(t)
S50 = fourier_series.truncate(50)
max_error = max(abs(f - S50).subs(t, t_val) for t_val in t_values)
print(f'Максимальная погрешность приближения: {max_error.evalf()}')
