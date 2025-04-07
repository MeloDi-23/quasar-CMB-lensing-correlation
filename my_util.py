"""
This util is a copy from corrfunc.util.
In the Corrfunc module, weighting is calculated using w1 w2; and if you specify weighting, an avg_weight will be returned.
However, when converting the count into wp and xi, the function does not make use of weighting at all!
I rewrite the convert_3d_counts_to_cf to make it make use of weight.
"""

def convert_3d_counts_to_cf(ND1, ND2, NR1, NR2,
                            D1D2, D1R2, D2R1, R1R2,
                            estimator='LS'):
    """
    Converts raw pair counts to a correlation function.

    Parameters
    ----------

    ND1 : integer | float
       Number of points in the first dataset
       Or sum of weights in the first dataset(if no weighting, this is the number of points.)
    ND2 : integer
        Number of points in the second dataset

    NR1 : integer
        Number of points in the randoms for first dataset

    NR2 : integer
        Number of points in the randoms for second dataset

    D1D2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and D2

    D1R2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and R2

    D2R1 : array-like, integer
        Pair-counts for the cross-correlation between D2 and R1

    R1R2 : array-like, integer
        Pair-counts for the cross-correlation between R1 and R2

    For all of these pair-counts arrays, the corresponding ``numpy``
    struct returned by the theory/mocks modules can also be passed

    estimator: string, default='LS' (Landy-Szalay)
        The kind of estimator to use for computing the correlation
        function. Currently, only supports Landy-Szalay

    Returns
    ---------

    cf : A numpy array
        The correlation function, calculated using the chosen estimator,
        is returned. NAN is returned for the bins where the ``RR`` count
        is 0.


    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from Corrfunc.theory.DD import DD
    >>> from Corrfunc.io import read_catalog
    >>> from Corrfunc.utils import convert_3d_counts_to_cf
    >>> X, Y, Z = read_catalog()
    >>> N = len(X)
    >>> boxsize = 420.0
    >>> rand_N = 3*N
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> rand_X = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Y = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Z = np.random.uniform(0, boxsize, rand_N)
    >>> nthreads = 2
    >>> rmin = 0.1
    >>> rmax = 15.0
    >>> nbins = 10
    >>> bins = np.linspace(rmin, rmax, nbins + 1)
    >>> autocorr = 1
    >>> DD_counts = DD(autocorr, nthreads, bins, X, Y, Z, boxsize=boxsize)
    >>> autocorr = 0
    >>> DR_counts = DD(autocorr, nthreads, bins,
    ...                X, Y, Z,
    ...                X2=rand_X, Y2=rand_Y, Z2=rand_Z, boxsize=boxsize)
    >>> autocorr = 1
    >>> RR_counts = DD(autocorr, nthreads, bins, rand_X, rand_Y, rand_Z,
    ...                boxsize=boxsize)
    >>> cf = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
    ...                              DD_counts, DR_counts,
    ...                              DR_counts, RR_counts)
    >>> for xi in cf: print("{0:10.6f}".format(xi))
    ...                    # doctest: +NORMALIZE_WHITESPACE
     22.769060
      3.612701
      1.621368
      1.000967
      0.691637
      0.511813
      0.398869
      0.318813
      0.255639
      0.207754

    """

    import numpy as np
    pair_counts = dict()
    weights = dict()
    fields = ['D1D2', 'D1R2', 'D2R1', 'R1R2']
    arrays = [D1D2, D1R2, D2R1, R1R2]
    for (field, array) in zip(fields, arrays):
        try:
            npairs = array['npairs']
            pair_counts[field] = npairs
        except IndexError:
            pair_counts[field] = array

        try:
            weight = array['weightavg']
            if np.isclose(weight, 0.0).all():
                # if no weight is provided, it will likely to return an weightavg of all 0.0s.
                weights[field] = np.ones_like(pair_counts[field])
            else:
                weights[field] = weight
        except IndexError:
            weights[field] = np.ones_like(pair_counts[field])

    pair_weight_product = {
        k: pair_counts[k]*weights[k]
        for k in fields
    }

    nbins = len(pair_counts['D1D2'])
    if (nbins != len(pair_counts['D1R2'])) or \
       (nbins != len(pair_counts['D2R1'])) or \
       (nbins != len(pair_counts['R1R2'])):
        msg = 'Pair counts must have the same number of elements (same bins)'
        raise ValueError(msg)

    nonzero = pair_counts['R1R2'] > 0
    if 'LS' in estimator or 'Landy' in estimator:
        fN1 = np.float64(NR1) / np.float64(ND1)
        fN2 = np.float64(NR2) / np.float64(ND2)
        cf = np.zeros(nbins)
        cf[:] = np.nan
        cf[nonzero] = (fN1 * fN2 * pair_weight_product['D1D2'][nonzero] -
                       fN1 * pair_weight_product['D1R2'][nonzero] -
                       fN2 * pair_weight_product['D2R1'][nonzero] +
                       pair_weight_product['R1R2'][nonzero]) / pair_weight_product['R1R2'][nonzero]
        if len(cf) != nbins:
            msg = 'Bug in code. Calculated correlation function does not '\
                  'have the same number of bins as input arrays. Input bins '\
                  '={0} bins in (wrong) calculated correlation = {1}'.format(
                      nbins, len(cf))
            raise RuntimeError(msg)
    else:
        msg = "Only the Landy-Szalay estimator is supported. Pass estimator"\
              "='LS'. (Got estimator = {0})".format(estimator)
        raise ValueError(msg)

    return cf


def convert_rp_pi_counts_to_wp(ND1, ND2, NR1, NR2,
                               D1D2, D1R2, D2R1, R1R2,
                               nrpbins, pimax, dpi=1.0,
                               estimator='LS'):
    """
    Converts raw pair counts to a correlation function.

    Parameters
    ----------

    ND1 : integer
       Number of points in the first dataset

    ND2 : integer
        Number of points in the second dataset

    NR1 : integer
        Number of points in the randoms for first dataset

    NR2 : integer
        Number of points in the randoms for second dataset

    D1D2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and D2

    D1R2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and R2

    D2R1 : array-like, integer
        Pair-counts for the cross-correlation between D2 and R1

    R1R2 : array-like, integer
        Pair-counts for the cross-correlation between R1 and R2

    For all of these pair-counts arrays, the corresponding ``numpy``
    struct returned by the theory/mocks modules can also be passed

    nrpbins : integer
        Number of bins in ``rp``

    pimax : float
        Integration distance along the line of sight direction

    dpi : float, default=1.0 Mpc/h
        Binsize in the line of sight direction

    estimator: string, default='LS' (Landy-Szalay)
        The kind of estimator to use for computing the correlation
        function. Currently, only supports Landy-Szalay

    Returns
    ---------

    wp : A numpy array
        The projected correlation function, calculated using the chosen
        estimator, is returned. If *any* of the ``pi`` bins (in an ``rp``
        bin) contains 0 for the ``RR`` counts, then ``NAN`` is returned
        for that ``rp`` bin.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from Corrfunc.theory.DDrppi import DDrppi
    >>> from Corrfunc.io import read_catalog
    >>> from Corrfunc.utils import convert_rp_pi_counts_to_wp
    >>> X, Y, Z = read_catalog()
    >>> N = len(X)
    >>> boxsize = 420.0
    >>> rand_N = 3*N
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> rand_X = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Y = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Z = np.random.uniform(0, boxsize, rand_N)
    >>> nthreads = 4
    >>> pimax = 40.0
    >>> nrpbins = 20
    >>> rpmin = 0.1
    >>> rpmax = 10.0
    >>> bins = np.linspace(rpmin, rpmax, nrpbins + 1)
    >>> autocorr = 1
    >>> DD_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    X, Y, Z, boxsize=boxsize)
    >>> autocorr = 0
    >>> DR_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    X, Y, Z,
    ...                    X2=rand_X, Y2=rand_Y, Z2=rand_Z, boxsize=boxsize)
    >>> autocorr = 1
    >>> RR_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    rand_X, rand_Y, rand_Z, boxsize=boxsize)
    >>> wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
    ...                                 DD_counts, DR_counts,
    ...                                 DR_counts, RR_counts,
    ...                                 nrpbins, pimax)
    >>> for w in wp: print("{0:10.6f}".format(w))
    ...                    # doctest: +NORMALIZE_WHITESPACE
    187.591897
     83.059026
     53.200243
     40.389026
     33.355778
     29.044893
     26.087995
     23.627759
     21.703655
     20.152961
     18.724304
     17.432795
     16.286740
     15.443105
     14.435802
     13.592479
     12.920796
     12.329687
     11.696258
     11.208016

    """

    import numpy as np
    if dpi <= 0.0:
        msg = 'Binsize along the line of sight (dpi) = {0}'\
              'must be positive'.format(dpi)
        raise ValueError(msg)

    xirppi = convert_3d_counts_to_cf(ND1, ND2, NR1, NR2,
                                     D1D2, D1R2, D2R1, R1R2,
                                     estimator=estimator)
    wp = np.empty(nrpbins)
    npibins = len(xirppi) // nrpbins
    if ((npibins * nrpbins) != len(xirppi)):
        msg = 'Number of pi bins could not be calculated correctly.'\
              'Expected to find that the total number of bins = {0} '\
              'would be the product of the number of pi bins = {1} '\
              'and the number of rp bins = {2}'.format(len(xirppi),
                                                       npibins,
                                                       nrpbins)
        raise ValueError(msg)

    # Check that dpi/pimax/npibins are consistent
    # Preventing issue #96 (https://github.com/manodeep/Corrfunc/issues/96)
    # where npibins would be calculated incorrectly, and the summation would
    # be wrong.
    if (dpi*npibins != pimax):
        msg = 'Pimax = {0} should be equal to the product of '\
              'npibins = {1} and dpi = {2}. Check your binning scheme.'\
              .format(pimax, npibins, dpi)
        raise ValueError(msg)

    for i in range(nrpbins):
        wp[i] = 2.0 * dpi * np.sum(xirppi[i * npibins:(i + 1) * npibins])

    return wp

