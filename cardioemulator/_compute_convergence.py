import numpy as np
import pandas as pd

def compute_convergence(file_PV_loops, period):
    """
    Evaluate the convergence to the limit cycle of a series of PV loops.

    Parameters
    ----------
    file_PV_loops : str
        Path of a csv file containing pressure and volume transients.
        Expected columns:

        - ``time``: time
        - ``VLV``: volume
        - ``pLV``: pressure

    period : float
        Heartbeat period (measure unit must be the same of the `time` column).

    Returns
    -------
    res : dict
        Dictionary with convergence metrics (listed below).
        Each metric contains a list, associated with the differences between
        consecutive heartbeats, or between a given heartbeat and the last one.
        If `n` PV loops are provided, each list contains `n - 1` elements.

        - ``err_rel_2_prs``: relative mean square error (norm 2) on pressure (between cycles)
        - ``err_rel_2_vol``: relative mean square error (norm 2) on volume (between cycles)
        - ``err_rel_I_prs``: relative maximum error (infinity) on pressure (between cycles)
        - ``err_rel_I_vol``: relative maximum error (infinity) on volume (between cycles)
        - ``err_rel_2_prs_E``: relative mean square error (norm 2) on pressure (wrt last cycle)
        - ``err_rel_2_vol_E``: relative mean square error (norm 2) on volume (wrt last cycle)
        - ``err_rel_I_prs_E``: relative maximum error (infinity) on pressure (wrt last cycle)
        - ``err_rel_I_vol_E``: relative maximum error (infinity) on volume (wrt last cycle)
        - ``err_rel_2_tot``: relative mean square error (norm 2) on pressure and volume (between cycles)
        - ``err_rel_2_tot_E``: relative mean square error (norm 2) on pressure and volume (wrt last cycle)
        - ``n_loops``: number of heartbeats

    """
    data_loop = pd.read_csv(file_PV_loops)

    err_rel_2_prs = list()
    err_rel_2_vol = list()
    err_rel_I_prs = list()
    err_rel_I_vol = list()

    err_rel_2_prs_E = list()
    err_rel_2_vol_E = list()
    err_rel_I_prs_E = list()
    err_rel_I_vol_E = list()

    err_rel_2_tot = list()
    err_rel_2_tot_E = list()

    iters_per_loop = sum(data_loop.time < period)
    n_loops = int(np.round(data_loop.time.max() / period))
    data_loops = [data_loop[(i + 0) * iters_per_loop:(i + 1) * iters_per_loop] for i in range(n_loops)]

    volE = np.array(data_loops[-1].VLV)
    prsE = np.array(data_loops[-1].pLV)
    nrm_2_prs_E = np.linalg.norm(prsE, 2     )
    nrm_2_vol_E = np.linalg.norm(volE, 2     )
    nrm_I_prs_E = np.linalg.norm(prsE, np.Inf)
    nrm_I_vol_E = np.linalg.norm(volE, np.Inf)

    for i in range(n_loops - 1):
        data1 = data_loops[i]
        data2 = data_loops[i+1]

        prs1 = np.array(data1.pLV)
        prs2 = np.array(data2.pLV)
        vol1 = np.array(data1.VLV)
        vol2 = np.array(data2.VLV)
        err_2_prs = np.linalg.norm(prs1 - prs2, 2     )
        nrm_2_prs = np.linalg.norm(prs2       , 2     )
        err_2_vol = np.linalg.norm(vol1 - vol2, 2     )
        nrm_2_vol = np.linalg.norm(vol2       , 2     )
        err_I_prs = np.linalg.norm(prs1 - prs2, np.Inf)
        nrm_I_prs = np.linalg.norm(prs2       , np.Inf)
        err_I_vol = np.linalg.norm(vol1 - vol2, np.Inf)
        nrm_I_vol = np.linalg.norm(vol2       , np.Inf)
        err_rel_2_prs.append(err_2_prs / nrm_2_prs)
        err_rel_2_vol.append(err_2_vol / nrm_2_vol)
        err_rel_I_prs.append(err_I_prs / nrm_I_prs)
        err_rel_I_vol.append(err_I_vol / nrm_I_vol)
        err_2_prs_E = np.linalg.norm(prs1 - prsE, 2     )
        err_2_vol_E = np.linalg.norm(vol1 - volE, 2     )
        err_I_prs_E = np.linalg.norm(prs1 - prsE, np.Inf)
        err_I_vol_E = np.linalg.norm(vol1 - volE, np.Inf)
        err_rel_2_prs_E.append(err_2_prs_E / nrm_2_prs_E)
        err_rel_2_vol_E.append(err_2_vol_E / nrm_2_vol_E)
        err_rel_I_prs_E.append(err_I_prs_E / nrm_I_prs_E)
        err_rel_I_vol_E.append(err_I_vol_E / nrm_I_vol_E)

        err_rel_2_tot.append(  err_2_prs   / nrm_2_prs   + err_2_vol   / nrm_2_vol  )
        err_rel_2_tot_E.append(err_2_prs_E / nrm_2_prs_E + err_2_vol_E / nrm_2_vol_E)

    results = dict()
    results['n_loops'] = n_loops
    results['err_rel_2_prs'] = err_rel_2_prs
    results['err_rel_2_vol'] = err_rel_2_vol
    results['err_rel_I_prs'] = err_rel_I_prs
    results['err_rel_I_vol'] = err_rel_I_vol
    results['err_rel_2_prs_E'] = err_rel_2_prs_E
    results['err_rel_2_vol_E'] = err_rel_2_vol_E
    results['err_rel_I_prs_E'] = err_rel_I_prs_E
    results['err_rel_I_vol_E'] = err_rel_I_vol_E
    results['err_rel_2_tot'] = err_rel_2_tot
    results['err_rel_2_tot_E'] = err_rel_2_tot_E

    return results
