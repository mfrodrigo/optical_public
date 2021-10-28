"""
Filter Design.
"""
import numpy as np
from numpy import ndarray, pi, abs, log10
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft


def transfer_function(frequencies: ndarray) -> ndarray:
    """
    Calculate the H(f) for a 4 order
    butterworth filter and cutoff frequency equal to 1 GHz.
    Args:
        frequencies: (ndarray) frequency vector

    Returns:
        H: (ndarray) Response in a certain region of filter frequencies.
    """
    s = 1j * 2 * pi * frequencies
    H = 1.559e39 / (s ** 4 + 1.642e10 * s ** 3 + 1.348e20 * s ** 2 + 6.482e29 * s + 1.559e39)
    return H


if __name__ == "__main__":
    f = np.linspace(0.1, 2, 1000) * 1e9
    H = transfer_function(f)
    mod_H = 20 * log10(abs(H))
    fig = plt.figure()
    plt.plot(f, mod_H, label='Função de transferência')
    plt.xlabel('frequência (Hz)')
    plt.ylabel('20*log10(|H|)')
    plt.legend()
    plt.grid(True)
    fig.savefig("resposta em frequência filtro.png")

    # fs = 10 GHz, f1 = 2 GHz e f2 = 0.5 GHz
    fs = 10e9
    f1 = 2e9
    f2 = 0.5e9
    N = int(400)
    t = np.linspace(0, N / fs, N, False)  # 1 second
    sig = np.sin(2 * pi * f1 * t) + np.sin(2 * pi * f2 * t)
    fig = plt.figure()
    plt.plot(t, sig, label='Sinal')
    plt.xlabel('tempo (s)')
    plt.ylabel('sin(2*pi*f1*t) + sin(2*pi*f2*t)')
    plt.legend()
    plt.grid(True)
    fig.savefig("sinal.png")

    sig_f = fft(sig)
    f = fftfreq(N, 1 / fs)
    fig = plt.figure()
    plt.plot(f[0:N // 2], 2.0 / N * np.abs(sig_f[0:N // 2]))
    plt.grid()
    plt.xlabel('frequência (Hz)')
    plt.ylabel('sin(2*pi*f1*t) + sin(2*pi*f2*t)')
    plt.grid(True)
    fig.savefig("fft.png")

    sig_filter = ifft(sig_f * transfer_function(f))
    fig = plt.figure()
    plt.plot(t, np.abs(sig_filter), label='Resposta ao filtro')
    plt.xlabel('tempo (s)')
    plt.ylabel('sin(2*pi*f1*t) + sin(2*pi*f2*t)')
    plt.legend()
    plt.grid(True)
    fig.savefig("sinal_filtrado.png")

    sig_f = fft(sig_filter)
    f = fftfreq(N, 1 / fs)
    fig = plt.figure()
    plt.plot(f[0:N // 2], 2.0 / N * np.abs(sig_f[0:N // 2]))
    plt.grid()
    plt.xlabel('frequência (Hz)')
    plt.ylabel('sin(2*pi*f1*t) + sin(2*pi*f2*t)')
    plt.grid(True)
    fig.savefig("fft_sinal_filtrado.png")