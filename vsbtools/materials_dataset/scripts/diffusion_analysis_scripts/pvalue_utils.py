import math
from typing import Callable, Dict


def _norm_cdf(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _validate_pvalue_inputs(callables, targets, margins, alternative):
    missing_targets = [k for k in callables.keys() if k not in targets]
    missing_margins = [k for k in callables.keys() if k not in margins]
    if missing_targets:
        raise KeyError(f"Missing targets for callable keys: {missing_targets}")
    if missing_margins:
        raise KeyError(f"Missing margins for callable keys: {missing_margins}")
    if alternative not in ("greater", "less", "two-sided"):
        raise ValueError("alternative must be one of: 'greater', 'less', 'two-sided'")


def _make_good_counter(callables, targets, margins):
    def is_good(entry):
        for key, fn in callables.items():
            if abs(fn(entry) - targets[key]) > margins[key]:
                return False
        return True

    def ds_good_count(ds):
        return sum(1 for e in ds if is_good(e))

    return ds_good_count


def _dataset_counts(ds_non_guided, ds_guided_local, ds_good_count):
    n0 = len(ds_non_guided)
    n1 = len(ds_guided_local)
    if n0 == 0 or n1 == 0:
        raise ValueError("Both datasets must be non-empty")
    k0 = ds_good_count(ds_non_guided)
    k1 = ds_good_count(ds_guided_local)
    return k0, n0, k1, n1


def _two_proportion_z_pvalue_from_counts(k0, n0, k1, n1, alternative: str):
    p0_hat = k0 / n0
    p1_hat = k1 / n1
    z_num = p1_hat - p0_hat
    p_pool = (k0 + k1) / (n0 + n1)
    se2 = p_pool * (1.0 - p_pool) * (1.0 / n0 + 1.0 / n1)
    if se2 <= 0.0:
        if z_num > 0.0:
            z_score = float("inf")
        elif z_num < 0.0:
            z_score = float("-inf")
        else:
            z_score = 0.0
        if alternative == "greater":
            p_value = 0.0 if z_num > 0.0 else 1.0
        elif alternative == "less":
            p_value = 0.0 if z_num < 0.0 else 1.0
        else:
            p_value = 0.0 if z_num != 0.0 else 1.0
        return z_score, float(p_value)

    z_score = z_num / math.sqrt(se2)
    if alternative == "greater":
        p_value = 1.0 - _norm_cdf(z_score)
    elif alternative == "less":
        p_value = _norm_cdf(z_score)
    else:
        p_value = 2.0 * min(_norm_cdf(z_score), 1.0 - _norm_cdf(z_score))
    return float(z_score), float(min(1.0, max(0.0, p_value)))


def get_two_proportion_z_test(
    ds_non_guided,
    callables: Dict[str, Callable],
    targets: Dict[str, float],
    margins: Dict[str, float],
    ds_guided=None,
    alternative: str = "greater",
):
    """
    Two-proportion pooled z-test on "good" shares.

    Returns:
        {"z_score": float, "p_value": float}
    """
    _validate_pvalue_inputs(callables, targets, margins, alternative)
    ds_good_count = _make_good_counter(callables, targets, margins)

    def _score_for(ds_guided_local):
        k0, n0, k1, n1 = _dataset_counts(ds_non_guided, ds_guided_local, ds_good_count)
        z_score, p_value = _two_proportion_z_pvalue_from_counts(k0, n0, k1, n1, alternative=alternative)
        return {"z_score": z_score, "p_value": p_value}

    if ds_guided is None:
        return _score_for
    return _score_for(ds_guided)


def get_p_value(
    ds_non_guided,
    callables: Dict[str, Callable],
    targets: Dict[str, float],
    margins: Dict[str, float],
    ds_guided=None,
    alternative: str = "greater",
    method: str = "exact",
    continuity_correction: bool = True,
):
    """
    Estimate guidance effect from the share of "good" entries.

    Entry is "good" if for every callable key:
        abs(callable(entry) - targets[key]) <= margins[key]

    If ds_guided is provided, returns a p-value (float) for the null hypothesis
    that guided and non-guided shares are equal.

    method:
        - "exact": exact one-sample binomial tail(s), with non-guided share as baseline
        - "normal": two-sample pooled-proportion z-test approximation

    If ds_guided is None, returns a scorer function that takes ds_guided and
    returns the corresponding p-value.
    """
    _validate_pvalue_inputs(callables, targets, margins, alternative)
    if method not in ("exact", "normal"):
        raise ValueError("method must be one of: 'exact', 'normal'")

    ds_good_count = _make_good_counter(callables, targets, margins)

    def _logsumexp(log_values):
        if not log_values:
            return float("-inf")
        max_log = max(log_values)
        if max_log == float("-inf"):
            return float("-inf")
        acc = sum(math.exp(v - max_log) for v in log_values)
        return max_log + math.log(acc)

    def _log_binom_pmf(k, n, p):
        if k < 0 or k > n:
            return float("-inf")
        if p <= 0.0:
            return 0.0 if k == 0 else float("-inf")
        if p >= 1.0:
            return 0.0 if k == n else float("-inf")
        # log(C(n,k) * p^k * (1-p)^(n-k)) using lgamma to avoid overflow.
        return (
            math.lgamma(n + 1)
            - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)
            + k * math.log(p)
            + (n - k) * math.log1p(-p)
        )

    def _one_sided_upper_tail(k, n, p):
        k = max(0, min(k, n + 1))
        if k >= n + 1:
            return 0.0
        log_tail = _logsumexp([_log_binom_pmf(i, n, p) for i in range(k, n + 1)])
        return float(min(1.0, max(0.0, math.exp(log_tail))))

    def _one_sided_lower_tail(k, n, p):
        k = max(-1, min(k, n))
        if k < 0:
            return 0.0
        log_tail = _logsumexp([_log_binom_pmf(i, n, p) for i in range(0, k + 1)])
        return float(min(1.0, max(0.0, math.exp(log_tail))))

    def _normal_approx_pooled(k0, n0, k1, n1):
        # Kept for API compatibility; no correction applied for pooled z-test.
        _ = continuity_correction
        _, p_value = _two_proportion_z_pvalue_from_counts(k0, n0, k1, n1, alternative=alternative)
        return p_value

    def _pvalue_for(ds_guided_local):
        k0, n0, k1, n1 = _dataset_counts(ds_non_guided, ds_guided_local, ds_good_count)
        p0 = k0 / n0

        if method == "normal":
            return _normal_approx_pooled(k0, n0, k1, n1)

        if alternative == "greater":
            return _one_sided_upper_tail(k1, n1, p0)
        if alternative == "less":
            return _one_sided_lower_tail(k1, n1, p0)
        p_low = _one_sided_lower_tail(k1, n1, p0)
        p_high = _one_sided_upper_tail(k1, n1, p0)
        return min(1.0, 2.0 * min(p_low, p_high))

    if ds_guided is None:
        return _pvalue_for
    return _pvalue_for(ds_guided)
