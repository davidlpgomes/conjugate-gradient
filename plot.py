import matplotlib.pyplot as plt
import math
import pandas as pd
import os


os.system('mkdir -p figures')

COLS = [
    'TIME',
    'L2_MISS_RATIO',
    'L3_BANDWIDTH',
    'FLOPS_DP',
    'FLOPS_AVX'
]


def plot(column, csv1, csv_old1, csv2, csv_old2):
    fig, ax = plt.subplots(figsize=(8, 5))

    abscissas = csv1['N'].values

    plt.plot(
        abscissas,
        csv1[column].values,
        label='OP1 - V2',
        linewidth=2,
        color='blue',
        antialiased=True
    )

    plt.plot(
        abscissas,
        csv_old1[column].values,
        label='OP1 - V1',
        linewidth=2,
        linestyle='dashed',
        color='red',
        antialiased=True
    )

    plt.plot(
        abscissas,
        csv2[column].values,
        label='OP2 - V2',
        linewidth=2,
        linestyle='dashdot',
        color='orange',
        antialiased=True
    )

    plt.plot(
        abscissas,
        csv_old2[column].values,
        label='OP2 - V1',
        linewidth=2,
        linestyle='dotted',
        color='green',
        antialiased=True
    )

    plt.legend(
        ncol=5, loc='best', fancybox=False, framealpha=1,
        borderpad=0.4, edgecolor='black'
    )

    plt.xscale('log')

    ax.set_xlabel('Tamanho do sistema linear')
    ax.set_ylabel(column)

    plt.title(f'{column}')
    plt.tight_layout()

    plt.savefig(f'figures/{column}.png', dpi=300)


def main():
    op1_wo = pd.read_csv('without_optimizations/op1.csv')
    op2_wo = pd.read_csv('without_optimizations/op2.csv')

    op1 = pd.read_csv('op1.csv')
    op2 = pd.read_csv('op2.csv')

    for c in COLS:
        plot(c, op1, op1_wo, op2, op2_wo)

    return


if __name__ == '__main__':
    main()
