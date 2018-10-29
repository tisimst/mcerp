import pytest
from ..lhd import ss, lhd


def test():
    # test single distribution
    d0 = ss.uniform(loc=-1, scale=2)  # uniform distribution,low=-1, width=2
    print(lhd(dist=d0, size=5))

    # test single distribution for multiple variables
    d1 = ss.norm(loc=0, scale=1)  # normal distribution, mean=0, stdev=1
    print(lhd(dist=d1, size=7, dims=5))

    # test multiple distributions
    d2 = ss.beta(2, 5)  # beta distribution, alpha=2, beta=5
    d3 = ss.expon(scale=1 / 1.5)  # exponential distribution, lambda=1.5
    print(lhd(dist=(d1, d2, d3), size=6))

    rand_lhs = lhd(dist=(d0, d1, d2, d3), size=100)
    spac_lhs = lhd(
        dist=(d0, d1, d2, d3),
        size=100,
        form="spacefilling",
        iterations=100,
        showcorrelations=True,
    )

    pytest.importorskip("matplotlib")
    pytest.importorskip("scatterplot_matrix")
    try:
        from scatterplot_matrix import scatterplot_matrix as spm
        import matplotlib.pyplot as plt
    except ImportError:
        print(rand_lhs)
        print(spac_lhs)
    else:
        names = ["U(-1,1)", "N(0,1)", "Beta(2,5)", "Exp(1.5)"]
        spm(rand_lhs.T, names=names)
        plt.suptitle("Completely Random LHS Design")
        plt.show()
        spm(spac_lhs.T, names=names)
        plt.suptitle("Space-Filling LHS Design")
        plt.show()
