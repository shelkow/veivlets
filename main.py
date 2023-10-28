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
Задание 2 (преобразование Фурье)
import math
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
gr = 5
num = 21
# нахождение сигнала f(t) = s(t) + p(t) + r(t) + n(t)
t = np.linspace(0, 2 * np.pi, 256)
s = t * (num + t) ** (1 / 2) # тренда s(t)
p = np.where(t < np.pi, np.cos(2 * (gr * t + num)), np.cos((gr * t + num))) #
периодические компоненты p(t)
r = np.random.uniform(-0.5, 0.5, 256) # случайные состовляющие r(t)
n = np.zeros(256) # импульсный шум n(t), равный нулю в остальных точках кроме
t=64, в это
n[64] = -8 # импульсный шум, в точке t=64 равен -8
f = s + p + r + n
fNotRand = s + p + n
# нахождение ряда фурье
X = np.fft.fft(f)
# построение графика амплитудного спектра
N = len(X) # длина сигнала
n = np.arange(N)
T = N / 256
freq = n / T # диапазон частот
plt.stem(freq, np.abs(X) / (N // 2), 'b', markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(0, 256)
plt.show()
# 20% наибольших по абсолютной величине коэффициентов
indexes = np.argsort(np.abs(X))
Abslen = int(len(X) * 0.2)
# Восстановление сигнала по оставшимся величинам c помощью обратного
преобразования фурье
newX = np.copy(X)
newX[indexes[:-Abslen]] = 0
reF = np.fft.ifft(newX)
# Разница между сигналами с помощью L1 нормы
print(f'\nРазница сигналов по L1 = {np.linalg.norm(fNotRand - reF, ord=1)}\n')
# Построение графиков исходного и восстановленного сигналов
plt.plot(t, fNotRand, 'g', )
plt.ylabel('Amplitude')
plt.show()
plt.plot(t, reF, 'r')
plt.ylabel('Amplitude')
plt.show()
print('Отображение на одном графике')
plt.plot(t, fNotRand, 'g', )
plt.plot(t, reF, 'r')
plt.ylabel('Amplitude')
plt.show()
Задание 3 (сжатие сигнала)
import math
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import pywt
from skimage.restoration import (denoise_wavelet, estimate_sigma)
gr = 5
num = 21
# нахождение сигнала f(t) = s(t) + p(t) + r(t) + n(t)
t = np.linspace(0, 2 * np.pi, 256)
s = t * (num + t) ** (1 / 2) # тренда s(t)
p = np.where(t < np.pi, 2 * np.cos(2 * (gr * t + num)), 3 * np.cos((gr * t +
num))) # периодические компоненты p(t)
r = np.random.uniform(-0.5, 0.5, 256) # случайные состовляющие r(t)
n = np.zeros(256) # импульсный шум n(t), равный нулю в остальных точках кроме
t=64, в это
n[64] = -8 # импульсный шум, в точке t=64 равен -8
f = s + p + r + n
fNotRand = s + p + n
# Вейвлет-преобразование Добеши2 print("Вейвлет-преобразование Добеши2:")
X = pywt.wavedec(f, 'db2', level=3)
theshhold = np.percentile(np.abs(np.concatenate(X[1:])), 80)
for i in range(1, len(X)): X[i] = pywt.threshold(X[i], theshhold, mode='soft')
refF = pywt.waverec(X, 'db2')
print(f'Разница сигналов по L1 = {np.linalg.norm(fNotRand - refF, ord=1)}')
plt.plot(t, fNotRand, 'b')
plt.plot(t, refF, 'r')
plt.ylabel('Добеши2')
plt.show()
# Вейвлет-преобразование Добеши10 print("Вейвлет-преобразование Добеши10:")
X = pywt.wavedec(f, 'db10', level=3)
theshhold = np.percentile(np.abs(np.concatenate(X[1:])), 80)
for i in range(1, len(X)): X[i] = pywt.threshold(X[i], theshhold, mode='soft')
refF = pywt.waverec(X, 'db10')
print(f'Разница сигналов по L1 = {np.linalg.norm(fNotRand - refF, ord=1)}')
plt.plot(t, fNotRand, 'b')
plt.plot(t, refF, 'r')
plt.ylabel('Добеши10')
plt.show()
# Вейвлет-преобразование биортагонального вейвлета 2,2 print("Вейвлетпреобразование биортагонального вейвлета 2,2:") X = pywt.wavedec(f, 'bior2.2',
level = 3)
X = pywt.wavedec(f, 'bior2.2', level=3)
theshhold = np.percentile(np.abs(np.concatenate(X[1:])), 80)
for i in range(1, len(X)): X[i] = pywt.threshold(X[i], theshhold, mode='soft')
refF = pywt.waverec(X, 'bior2.2')
print(f'Разница сигналов по L1 = {np.linalg.norm(fNotRand - refF, ord=1)}')
plt.plot(t, fNotRand, 'b')
plt.plot(t, refF, 'r')
plt.ylabel('Биортогональный вейвлет 2,2')
plt.show()
# Вейвлет-преобразование симметричного вейвлета print("Вейвлет-преобразование
симметричного вейвлета:")
X = pywt.wavedec(f, 'sym4', level=3)
theshhold = np.percentile(np.abs(np.concatenate(X[1:])), 80)
for i in range(1, len(X)): X[i] = pywt.threshold(X[i], theshhold, mode='soft')
refF = pywt.waverec(X, 'sym4')
print(f'Разница сигналов по L1 = {np.linalg.norm(fNotRand - refF, ord=1)}')
plt.plot(t, fNotRand, 'b')
plt.plot(t, refF, 'r')
plt.ylabel('Симметричный вейвлет 4')
plt.show()
Задание 4 (устранение шума в сигнале)
import math
import numpy as np
import matplotlib.pyplot as plt
import pywt
from skimage.restoration import (denoise_wavelet, estimate_sigma)
gr = 5
num = 21
t = np.linspace(0, 2 * np.pi, 256)
s = t * (num + t) ** (1 / 2)
p = np.where(t < np.pi, 2 * np.cos(2 * (gr * t + num)), 3 * np.cos((gr * t +
num)))
n = np.random.normal(0, 4, 256)
fNoise = s + p + n
f = s + p
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
print(-20 * np.log10(np.linalg.norm(abs(fNoise - f)) / np.linalg.norm(f)))
fDenoise = denoise_wavelet(fNoise, wavelet='db2')
print(f"значение соотношения без шума", - 20 *
np.log10(np.linalg.norm(abs(fDenoise - f)) / np.linalg.norm(f)))
axs.flat[0].plot(t, f, 'b', label="Исходный сигнал")
axs.flat[0].plot(t, fNoise, 'g', label="Исходный сигнал (с шумом)")
axs.flat[0].plot(t, fDenoise, 'r', label="Восстановленный сигнал")
axs.flat[0].legend()
axs.flat[0].set_title('Добеши2')
fDenoise = denoise_wavelet(fNoise, wavelet='db10')
print(f"значение соотношения без шума:", - 20 *
np.log10(np.linalg.norm(abs(fDenoise - f)) / np.linalg.norm(f)))
axs.flat[1].plot(t, f, 'b', label="Исходный сигнал")
axs.flat[1].plot(t, fNoise, 'g', label="Исходный сигнал (с шумом)")
axs.flat[1].plot(t, fDenoise, 'r', label="Восстановленный сигнал")
axs.flat[1].legend()
axs.flat[1].set_title('Добеши10')
fDenoise = denoise_wavelet(fNoise, wavelet='bior2.2')
print(f"sзначение соотношения без шума", - 20 *
np.log10(np.linalg.norm(abs(fDenoise - f)) / np.linalg.norm(f)))
axs.flat[2].plot(t, f, 'b', label="Исходный сигнал")
axs.flat[2].plot(t, fNoise, 'g', label="Исходный сигнал (с шумом)")
axs.flat[2].plot(t, fDenoise, 'r', label="Восстановленный сигнал")
axs.flat[2].legend()
axs.flat[2].set_title('Биортогональный вейвлет 2,2')
fDenoise = denoise_wavelet(fNoise, wavelet='sym4')
print(f"значение соотношения без шума", - 20 *
np.log10(np.linalg.norm(abs(fDenoise - f)) / np.linalg.norm(f)))
axs.flat[3].plot(t, f, 'b', label="Исходный сигнал")
axs.flat[3].plot(t, fNoise, 'g', label="Исходный сигнал (с шумом)")
axs.flat[3].plot(t, fDenoise, 'r', label="Восстановленный сигнал")
axs.flat[3].legend()
axs.flat[3].set_title('Симметричный вейвлет 4')
