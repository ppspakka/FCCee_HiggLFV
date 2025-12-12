import os
import glob
from typing import List, Dict, Any, Union, Optional

import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.colors import LogNorm

# Numpy-backed histogram wrappers (project-specific)
from numpy_hist import NumpyHist1D, NumpyHist2D

__all__ = [
    "find_histogram_names",
    "find_histogram_names_multiple_keywords",
    "read_root_histograms",
    "rebinning_histogram",
    "apply_weights",
    "plot_target_histogram",
    "merge_histograms_by_mapping",
    "plot_2d_histogram",
    "plot_brazil_limits",
    "upper_limit_xsec_br",
    "optimize_cut_punzi",
    "feldman_cousins_upper_limit",
]


def find_histogram_names(all_names: List[str],
                         keyword: str,
                         only_one: bool = True,
                         mode: str = "contains",
                         ignore_keyword: Optional[str] = None,
                         get_list: bool = False) -> Union[str, List[str]]:
    """Find histogram name(s) in `all_names` matching `keyword`.

    mode: 'contains' or 'equals'.
    If get_list is True returns list of matches.
    If only_one is True and not exactly one match is found, raises ValueError.
    """
    if mode == "contains":
        matched = [n for n in all_names if keyword in n]
    elif mode == "equals":
        matched = [n for n in all_names if keyword == n]
    else:
        raise ValueError("mode must be 'contains' or 'equals'")

    if ignore_keyword is not None:
        matched = [n for n in matched if ignore_keyword not in n]

    if get_list:
        return matched

    if only_one:
        if len(matched) == 1:
            return matched[0]
        raise ValueError(f"Expected exactly one match for keyword '{keyword}', found {len(matched)}.")
    return matched


def find_histogram_names_multiple_keywords(all_names: List[str],
                                           keywords: List[str],
                                           only_one: bool = True,
                                           mode: str = "contains") -> Union[str, List[str]]:
    """Find histogram names matching ALL of the provided keywords."""
    # matched = []
    # for kw in keywords:
    #     if mode == "contains":
    #         matched.extend([n for n in all_names if kw in n])
    #     elif mode == "equals":
    #         matched.extend([n for n in all_names if kw == n])
    #     else:
    #         raise ValueError("mode must be 'contains' or 'equals'")

    # # unique
    # matched = list(set(matched))

    # if only_one:
    #     if len(matched) == 1:
    #         return matched[0]
    #     raise ValueError(f"Expected exactly one match for keywords '{keywords}', found {len(matched)}.")
    # return matched
    matched = all_names.copy()
    for kw in keywords:
        if mode == "contains":
            matched = [n for n in matched if kw in n]
        elif mode == "equals":
            matched = [n for n in matched if kw == n]
        else:
            raise ValueError("mode must be 'contains' or 'equals'")

    if only_one:
        if len(matched) == 1:
            return matched[0]
        raise ValueError(f"Expected exactly one match for keywords '{keywords}', found {len(matched)}.")
    return matched


def read_root_histograms(outdir: str = "outdir",
                         pattern: str = "*.root",
                         blacklist: List[str] = []) -> Dict[str, Dict[str, Any]]:
    """Read top-level histograms from ROOT files and convert to numpy-backed wrappers.

    Returns dict: { file_base_name: { histname: NumpyHist1D/2D, ... }, ... }
    """
    result: Dict[str, Dict[str, Any]] = {}
    ROOT.gROOT.SetBatch(True)

    paths = sorted(glob.glob(os.path.join(outdir, pattern)))
    print(f"Found {len(paths)} ROOT files matching pattern '{pattern}' in '{outdir}'.")

    for path in paths:
        print(path, flush=True)
        if any(bl in path for bl in blacklist):
            continue

        fname = os.path.splitext(os.path.basename(path))[0]
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            if f:
                f.Close()
            continue

        result[fname] = {}
        keys = f.GetListOfKeys()
        if not keys:
            f.Close()
            continue

        for key in keys:
            obj = key.ReadObj()
            if not obj:
                continue

            # check TH2 before TH1 because TH2 inherits from TH1
            if obj.InheritsFrom("TH2"):
                try:
                    h = NumpyHist2D.from_root_th2(obj)
                except Exception as e:
                    print(f"Warning: failed to convert 2D histogram '{obj.GetName()}': {e}")
                    continue
                result[fname][obj.GetName()] = h
            elif obj.InheritsFrom("TH1"):
                try:
                    h = NumpyHist1D.from_root_th1(obj)
                except Exception as e:
                    print(f"Warning: failed to convert histogram '{obj.GetName()}': {e}")
                    continue
                result[fname][obj.GetName()] = h

        f.Close()
    return result


def rebinning_histogram(hist, rebin_factor: Union[int, tuple, list, np.ndarray]):
    """Rebin a histogram by `rebin_factor`. Returns a clone (new object)."""
    # If factor is iterable, check whether all entries <= 1 (no-op)
    if isinstance(rebin_factor, (list, tuple, np.ndarray)):
        rf_is_one = all(int(x) <= 1 for x in rebin_factor)
    else:
        rf_is_one = int(rebin_factor) <= 1

    if rf_is_one:
        return hist.Clone()

    rebinned_hist = hist.Clone()
    # preserve the original approach: call Rebin (may accept tuples in your wrappers)
    rebinned_hist.Rebin(rebin_factor)
    return rebinned_hist


def apply_weights(hist_dict: Dict[str, Dict[str, Any]],
                  cross_sections: Dict[str, float],
                  target_lumi: float,
                  all_hist_lists: Optional[List[str]] = None):
    """Apply per-process weight = xsec * lumi / N_generated (taken from initial histogram).

    If `all_hist_lists` is provided it is used by name lookups; otherwise caller's global
    `all_hist_lists` will be used if present.
    """
    for proc_name, histograms in hist_dict.items():
        if proc_name not in cross_sections:
            continue
        xsec = cross_sections[proc_name]

        # get initial histogram and the final selection histogram for efficiency
        initial_hist = histograms.get("00_Initial_n_muons", None)

        # try to use provided all_hist_lists, otherwise try to find in globals (backwards compatible)
        search_hist_list = all_hist_lists if all_hist_lists is not None else globals().get("all_hist_lists", [])
        try:
            final_hist_name = find_histogram_names(search_hist_list, "finalstate_nocut_m_collinear", only_one=True)
        except Exception:
            final_hist_name = None
        last_hist = histograms.get(final_hist_name, None) if final_hist_name else None

        if initial_hist is None:
            continue

        total_events = initial_hist.GetEntries()
        if total_events == 0:
            continue

        weight = (xsec * target_lumi) / total_events
        passed_events = last_hist.GetEntries() if last_hist else 0
        efficiency = passed_events / total_events if total_events > 0 else 0
        total_yield = passed_events * weight

        print(f"Process: {proc_name}, Weight: {weight}, Total Events: {total_events}, "
              f"XSec: {xsec} pb, Efficiency: {efficiency:.2%}, Total Yield: {total_yield}")

        # scale all histograms for this process
        for hname, hist in histograms.items():
            hist.Scale(weight)

def plot_target_histogram(signal, background, histname,
                          custom_parameters={},
                          outdir="prefit_histograms",
                          saveas_pdf=False,
                          save_name="histogram_plot",
                          lumi_label="1 ab$^{-1}$ ($\\sqrt{s} = 240$ GeV)",
                          process_colours=None,
                          process_propernames=None,
                          figsize=(10,6)):
    """
    Plot one target histogram `histname` given:
      - signal:    {'signal_file': {'histname': TH1/TH2-like, ...}, ...}
      - background:{'bkg_file':    {'histname': TH1/TH2-like, ...}, ...}
    Produces a stacked background and overlaid signal(s). Handles absent inputs.
    (This function expects 1D histograms for plotting; 2D histograms should be handled elsewhere.)
    """

    custom_title = custom_parameters.get('title', None)
    custom_xlabel = custom_parameters.get('xlabel', None)
    custom_ylabel = custom_parameters.get('ylabel', None)
    y_unit = custom_parameters.get('yunit', 'unit')
    yscale = custom_parameters.get('yscale', 'linear')
    xlim = custom_parameters.get('xlim', None)
    ylim = custom_parameters.get('ylim', None)
    if yscale not in ['linear', 'log']:
        raise ValueError(f"Invalid yscale '{yscale}'; must be 'linear' or 'log'.")

    # collect any histogram to extract binning
    def get_hist_from_dict(d):
        for proc in d.values():
            if histname in proc:
                return proc[histname]
        return None

    ref_hist = get_hist_from_dict(background) or get_hist_from_dict(signal)
    if ref_hist is None:
        raise ValueError(f"No histogram named '{histname}' found in signal or background inputs.")

    # Guard: this drawer expects 1D; skip 2D with a clear message
    if hasattr(ref_hist, 'GetNbinsY') and callable(ref_hist.GetNbinsY) and ref_hist.GetNbinsY() > 1 and ref_hist.GetNbinsX() > 1:
        raise ValueError("plot_target_histogram currently supports 1D histograms only. Please use a 2D-specific plotting routine.")

    # If none, get the title from histogram title
    if custom_title is None:
        custom_title = ref_hist.GetTitle()
    if custom_xlabel is None:
        custom_xlabel = ref_hist.GetXaxis().GetTitle()
    if custom_ylabel is None:
        custom_ylabel = ref_hist.GetYaxis().GetTitle()

    # get bin edges and number of bins
    xa = ref_hist.GetXaxis()
    nbins = xa.GetNbins()
    edges = [xa.GetBinLowEdge(i) for i in range(1, nbins+1)]
    edges.append(xa.GetBinUpEdge(nbins))
    edges = np.array(edges)
    
    # Single bin width (assume single equal width for all bins)
    bin_width = int(edges[1] - edges[0])

    # helper to convert hist -> numpy arrays (bin contents and bin errors)
    def hist_to_arrays(h):
        contents = np.array([h.GetBinContent(i) for i in range(1, nbins+1)])
        errors   = np.array([h.GetBinError(i)   for i in range(1, nbins+1)])
        return contents, errors

    # prepare background merged per-process (top-level keys)
    procs = sorted(background.keys()) if background else []
    if not procs:
        procs = []


    merged_bkg = {}
    merged_bkg_err = {}
    for p in procs:
        proc_hist = background.get(p, {}).get(histname, None)
        if proc_hist is None:
            merged_bkg[p] = np.zeros(nbins)
            merged_bkg_err[p] = np.zeros(nbins)
        else:
            c, e = hist_to_arrays(proc_hist)
            merged_bkg[p] = c
            merged_bkg_err[p] = e

    proc_counts = {p: merged_bkg[p] for p in procs}
    sorted_procs = sorted(proc_counts.keys(), key=lambda p: np.sum(proc_counts[p]), reverse=False)

    # default colours/names if not provided
    if process_colours is None:
        cmap = plt.get_cmap("tab20")
        process_colours = {p: cmap(i % 20) for i, p in enumerate(sorted_procs)}
    if process_propernames is None:
        process_propernames = {p: p for p in procs}
        
    # total background error (sum in quadrature)
    if procs:
        total_err = np.sqrt(np.sum([merged_bkg_err[p]**2 for p in procs], axis=0))
    else:
        total_err = np.zeros(nbins)

    # prepare signals (may be multiple top-level keys)
    sig_procs = sorted(signal.keys()) if signal else []
    signals = {}
    for s in sig_procs:
        sh = signal.get(s, {}).get(histname, None)
        if sh is None:
            signals[s] = (np.zeros(nbins), np.zeros(nbins))
        else:
            signals[s] = hist_to_arrays(sh)

    # start plotting
    plt.figure(figsize=figsize)
    baseline = np.zeros(nbins)
    
    for p in sorted_procs:
        counts = merged_bkg[p]
        color = process_colours.get(p, None)
        label = process_propernames.get(p, p)
        # stairs(values, edges, baseline=baseline) will fill between baseline and values
        plt.stairs(baseline + counts, edges, baseline=baseline, label=label,
                   color=color, fill=True)
        baseline = baseline + counts
    

    # draw background uncertainty band (hatched)
    if procs:
        upper = baseline + total_err
        lower = baseline - total_err
        plt.fill_between(edges, np.append(lower, lower[-1]), np.append(upper, upper[-1]),
                         step='post', facecolor='none', hatch='///', edgecolor='black', linewidth=0,
                         label="background prefit unc.")

    # draw signal(s) as line (no fill)
    if sig_procs:
        line_styles = ['-', '--', '-.', ':']
        for i, s in enumerate(sig_procs):
            counts, errs = signals[s]
            if np.all(counts == 0):
                continue
            linestyle = line_styles[i % len(line_styles)]
            plt.stairs(counts, edges, label=s, color='black', linewidth=2, linestyle=linestyle)

    # if no background but signals exist, ensure baseline 0 and limits auto
    plt.xlim(edges[0], edges[-1])
    plt.xlabel(custom_xlabel)
    plt.ylabel(custom_ylabel + f' / {bin_width} {y_unit}')
    plt.yscale(yscale)
    # plt.legend(loc='best', bbox_to_anchor=(0, 0.5, 1, 0.5))
    # If number of legend entries is large, use multiple columns
    n_legend_cols = 1
    n_entries = len(procs) + len(sig_procs) + (1 if procs else 0)  # +1 for uncertainty band
    if n_entries > 10:
        n_legend_cols = 2
    plt.legend(loc='best', ncol=n_legend_cols)
    
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    # CMS-style label
    try:
        title = custom_title if custom_title is not None else histname
        hep.cms.label(exp=None, llabel=title + " , Preliminary", rlabel=lumi_label)
    except Exception:
        pass

    if saveas_pdf:
        # Ensure output directory exists
        os.makedirs(outdir, exist_ok=True)
        savepath = os.path.join(outdir, f"{save_name}.pdf")
        # plt.savefig(savepath)
        # Use bbox_inches='tight' to avoid cutting off labels
        plt.savefig(savepath, bbox_inches='tight')
        print(f"Saved plot to {savepath}")

    # Show plot instead of saving
    plt.show()


def merge_histograms_by_mapping(histogram_dict, mapping_dict):
    merged_dict = {}
    for old_name, histograms in histogram_dict.items():
        if old_name in mapping_dict:
            group_name = mapping_dict[old_name]
            if group_name not in merged_dict:
                merged_dict[group_name] = {}
            for hname, hist in histograms.items():
                if hname not in merged_dict[group_name]:
                    merged_dict[group_name][hname] = hist.Clone()
                else:
                    merged_dict[group_name][hname] = merged_dict[group_name][hname].Merge(hist)
        else:
            # No mapping, keep original
            merged_dict[old_name] = histograms
    return merged_dict

def merge_histograms_by_mapping(histogram_dict: Dict[str, Dict[str, Any]],
                                mapping_dict: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    """Merge histograms by mapping old_name -> group_name.

    For collisions, use Merge method if available, otherwise overwrite with clone.
    """
    merged_dict: Dict[str, Dict[str, Any]] = {}
    for old_name, histograms in histogram_dict.items():
        if old_name in mapping_dict:
            group_name = mapping_dict[old_name]
            if group_name not in merged_dict:
                merged_dict[group_name] = {}
            for hname, hist in histograms.items():
                if hname not in merged_dict[group_name]:
                    merged_dict[group_name][hname] = hist.Clone()
                else:
                    # attempt Merge, fallback to adding contents if Merge returns None
                    merged = merged_dict[group_name][hname].Merge(hist)
                    if merged is None:
                        merged_dict[group_name][hname] = hist.Clone()
        else:
            merged_dict[old_name] = histograms
    return merged_dict

def plot_2d_histogram(hist,
                      custom_parameters={},
                      outdir="prefit_histograms",
                      saveas_pdf=False,
                      save_name="hist2d_plot",
                      lumi_label="1 ab$^{-1}$ ($\\sqrt{s} = 240$ GeV)",
                      cmap="viridis",
                      mask_zeros=True,
                      zscale="linear",
                      contour=True,
                      contour_levels=None,
                      figsize=(10, 6)):
    """
    Plot a single 2D histogram (ROOT.TH2-like or NumpyHist2D-like) in a CMS-style panel.

    Inputs:
      - hist: a single 2D histogram object (ROOT TH2 or the project's NumpyHist2D wrapper).
      - custom_parameters: dict with overrides: title, xlabel, ylabel, zlabel, xlim, ylim.
      - mask_zeros: replace zero/negative bins with nan so they are not shown.
      - zscale: 'linear' or 'log' (applied to color scale).
      - contour: whether to draw contour lines on top of the color map.
      - contour_levels: explicit contour levels (list) or None to auto-set.
    Returns:
      - matplotlib Axes
    """
    from matplotlib.colors import LogNorm
    import matplotlib.pyplot as _plt

    title = custom_parameters.get("title", hist.GetTitle() if hasattr(hist, "GetTitle") else "")
    xlabel = custom_parameters.get("xlabel", hist.GetXaxis().GetTitle() if hasattr(hist, "GetXaxis") else "")
    ylabel = custom_parameters.get("ylabel", hist.GetYaxis().GetTitle() if hasattr(hist, "GetYaxis") else "")
    zlabel = custom_parameters.get("zlabel", "Counts")
    xlim = custom_parameters.get("xlim", None)
    ylim = custom_parameters.get("ylim", None)

    # Helper: extract bin edges and contents for a ROOT.TH2-like interface
    def _extract_from_root_th2(h):
        nx = int(h.GetNbinsX())
        ny = int(h.GetNbinsY())
        x_edges = [h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 1)]
        x_edges.append(h.GetXaxis().GetBinUpEdge(nx))
        y_edges = [h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 1)]
        y_edges.append(h.GetYaxis().GetBinUpEdge(ny))
        arr = np.zeros((nx, ny), dtype=float)
        for ix in range(1, nx + 1):
            for iy in range(1, ny + 1):
                arr[ix - 1, iy - 1] = h.GetBinContent(ix, iy)
        return np.array(x_edges), np.array(y_edges), arr

    # Helper: extract from numpy-backed hist wrappers (common names)
    def _extract_from_numpy_wrapper(h):
        # try common attributes: counts / values and xedges / x_edges
        if hasattr(h, "counts") and hasattr(h, "xedges") and hasattr(h, "yedges"):
            return np.array(h.xedges), np.array(h.yedges), np.array(h.counts)
        if hasattr(h, "values") and hasattr(h, "edges"):
            # values shape likely (nx, ny), edges is list [xedges, yedges]
            edges = h.edges
            return np.array(edges[0]), np.array(edges[1]), np.array(h.values)
        # fallback: try ROOT-like API on wrapper
        if hasattr(h, "GetNbinsX"):
            return _extract_from_root_th2(h)
        raise ValueError("Unrecognized 2D histogram object: couldn't extract arrays.")

    # extract arrays
    try:
        if hasattr(hist, "GetNbinsX") and hasattr(hist, "GetNbinsY"):
            x_edges, y_edges, contents = _extract_from_root_th2(hist)
        else:
            x_edges, y_edges, contents = _extract_from_numpy_wrapper(hist)
    except Exception as e:
        raise RuntimeError(f"Failed to extract 2D histogram arrays: {e}")

    # contents is shape (nx, ny) with x index first; for plotting we will transpose
    # mask zeros / negative values if requested
    if mask_zeros:
        contents_masked = contents.copy()
        contents_masked[contents_masked <= 0] = np.nan
    else:
        contents_masked = contents.copy()

    # Prepare plotting
    fig, ax = _plt.subplots(figsize=figsize)

    # Choose normalization for color scale
    norm = None
    if zscale == "log":
        # avoid passing LogNorm if all values are nan/zero
        if np.nanmax(contents_masked) > 0:
            norm = LogNorm(vmin=max(1e-12, np.nanmin(contents_masked[np.isfinite(contents_masked)])),
                           vmax=np.nanmax(contents_masked))
        else:
            norm = None

    # pcolormesh expects arrays of shape (nx+1, ny+1) for edges and data shape (nx, ny) transposed for display
    pcm = ax.pcolormesh(x_edges, y_edges, contents_masked.T,
                        cmap=cmap, shading="auto", norm=norm)

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(zlabel)

    # Optionally draw contours (use centers for contour placement)
    if contour:
        # compute centers and safe array for contour (replace nan with 0)
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        Z = contents_masked.copy().T  # shape (ny, nx) -> after transpose (ny, nx) appropriate for contour with meshgrid

        # choose contour levels automatically if not provided
        if contour_levels is None:
            finite = Z[np.isfinite(Z)]
            if finite.size > 0:
                # use percentiles for 1/2 sigma like bands or simple multiples
                low = np.nanpercentile(finite, 70)
                high = np.nanpercentile(finite, 95)
                contour_levels = sorted(list({low, high, np.nanmax(finite)}))
            else:
                contour_levels = []

        if len(contour_levels):
            Xg, Yg = np.meshgrid(x_centers, y_centers)
            try:
                cs = ax.contour(Xg, Yg, Z, levels=contour_levels, colors="k", linewidths=0.8)
                ax.clabel(cs, fmt="%1.2g", fontsize=8)
            except Exception:
                # ignore contour failures (e.g. insufficient finite points)
                pass

    # Labels, limits, title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # if title:
    #     ax.set_title(title)

    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(float(x_edges[0]), float(x_edges[-1]))
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(float(y_edges[0]), float(y_edges[-1]))

    # CMS-style label (best-effort)
    try:
        hep.cms.label(exp=None, llabel=(title if title else "") + " , Preliminary", rlabel=lumi_label, ax=ax)
    except Exception:
        try:
            hep.cms.label(exp=None, llabel=(title if title else "") + " , Preliminary", rlabel=lumi_label)
        except Exception:
            pass

    # Save if requested
    if saveas_pdf:
        os.makedirs(outdir, exist_ok=True)
        savepath = os.path.join(outdir, f"{save_name}.pdf")
        fig.savefig(savepath, bbox_inches="tight")
        print(f"Saved 2D plot to {savepath}")

    _plt.show()
    return ax



def plot_brazil_limits(data: Any,
                       mass_width: float = 2.0,
                       ax=None,
                       title: Optional[str] = None,
                       xlabel: str = "Mass (GeV)",
                       ylabel: str = "Upper limit",
                       show: bool = True,
                       color_median: str = "black"):
    """Plot 'Brazil' style expected-limit bands. Accepts several input formats (see docstring in notebook)."""
    import numpy as _np
    import matplotlib.pyplot as _plt

    # normalize input
    if isinstance(data, dict) and "mass" in data and "limits" in data:
        masses = _np.atleast_1d(data["mass"]).astype(float)
        raw_limits = list(data["limits"])
    else:
        masses = _np.atleast_1d(125.0).astype(float)
        raw_limits = [data]

    def _standardize(entry):
        out = {"median": None, "up1": None, "down1": None, "up2": None, "down2": None}
        if entry is None:
            return out
        if isinstance(entry, (int, float, _np.floating, _np.integer)):
            out["median"] = float(entry)
            return out
        if isinstance(entry, (list, tuple, _np.ndarray)):
            arr = list(entry)
            if len(arr) == 1:
                out["median"] = float(arr[0])
            elif len(arr) == 3:
                out["median"] = float(arr[0]); out["down1"] = float(arr[1]); out["up1"] = float(arr[2])
            elif len(arr) >= 5:
                out["median"] = float(arr[0]); out["down1"] = float(arr[1]); out["up1"] = float(arr[2])
                out["down2"] = float(arr[3]); out["up2"] = float(arr[4])
            return out
        if isinstance(entry, dict):
            if "0.5" in entry:
                out["median"] = float(entry.get("0.5"))
                out["up1"] = float(entry.get("0.84")) if "0.84" in entry else None
                out["down1"] = float(entry.get("0.16")) if "0.16" in entry else None
                out["up2"] = float(entry.get("0.975")) if "0.975" in entry else None
                out["down2"] = float(entry.get("0.025")) if "0.025" in entry else None
                return out
            if "median" in entry:
                out["median"] = float(entry.get("median"))
                out["up1"] = float(entry.get("up1")) if entry.get("up1") is not None else entry.get("+1", None)
                out["down1"] = float(entry.get("down1")) if entry.get("down1") is not None else entry.get("-1", None)
                out["up2"] = float(entry.get("up2")) if entry.get("up2") is not None else entry.get("+2", None)
                out["down2"] = float(entry.get("down2")) if entry.get("down2") is not None else entry.get("-2", None)
                return out
            for k, v in entry.items():
                lk = str(k).lower()
                try:
                    val = float(v)
                except Exception:
                    continue
                if "med" in lk:
                    out["median"] = val
                if lk in ("up1", "+1", "p1"):
                    out["up1"] = val
                if lk in ("down1", "-1", "m1"):
                    out["down1"] = val
                if lk in ("up2", "+2", "p2"):
                    out["up2"] = val
                if lk in ("down2", "-2", "m2"):
                    out["down2"] = val
            return out
        return out

    std_limits = [_standardize(x) for x in raw_limits]

    # sort by mass
    order = _np.argsort(masses)
    masses = masses[order]
    std_limits = [std_limits[i] for i in order]

    xs = _np.array(masses)
    medians = _np.array([L["median"] if L["median"] is not None else _np.nan for L in std_limits], dtype=float)
    up1 = _np.array([L["up1"] if L["up1"] is not None else _np.nan for L in std_limits], dtype=float)
    down1 = _np.array([L["down1"] if L["down1"] is not None else _np.nan for L in std_limits], dtype=float)
    up2 = _np.array([L["up2"] if L["up2"] is not None else _np.nan for L in std_limits], dtype=float)
    down2 = _np.array([L["down2"] if L["down2"] is not None else _np.nan for L in std_limits], dtype=float)

    if ax is None:
        fig, ax = _plt.subplots(figsize=(10, 8))

    if xs.size == 1:
        m = xs[0]
        x_left = m - mass_width
        x_right = m + mass_width
        if not _np.isnan(up2) and not _np.isnan(down2):
            ax.fill_between([x_left, x_right], [down2[0], down2[0]], [up2[0], up2[0]],
                            color="yellow", edgecolor="none", label=r"$\pm2\sigma$ (expected)", alpha=0.9)
        if not _np.isnan(up1) and not _np.isnan(down1):
            ax.fill_between([x_left, x_right], [down1[0], down1[0]], [up1[0], up1[0]],
                            color="green", edgecolor="none", label=r"$\pm1\sigma$ (expected)", alpha=0.9)
        if not _np.isnan(medians[0]):
            ax.hlines(medians[0], x_left, x_right, colors=color_median, linewidth=2, label="median expected")
            ax.plot(m, medians[0], marker="o", color=color_median)
    else:
        if _np.any(~_np.isnan(up2)) and _np.any(~_np.isnan(down2)):
            ax.fill_between(xs, down2, up2, where=~(_np.isnan(down2) | _np.isnan(up2)),
                            color="yellow", edgecolor="none", label=r"$\pm2\sigma$ (expected)", alpha=0.9, interpolate=True)
        if _np.any(~_np.isnan(up1)) and _np.any(~_np.isnan(down1)):
            ax.fill_between(xs, down1, up1, where=~(_np.isnan(down1) | _np.isnan(up1)),
                            color="green", edgecolor="none", label=r"$\pm1\sigma$ (expected)", alpha=0.9, interpolate=True)
        if not _np.all(_np.isnan(medians)):
            ax.plot(xs, medians, color=color_median, linestyle="-", linewidth=1, label="median expected")
            ax.plot(xs, medians, "o", color=color_median)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    hep.cms.label(exp="", llabel="Fast Simulation", rlabel="1 ab$^{-1}$ ($\\sqrt{s} = 240$ GeV)")
    xmin = float(_np.min(xs) - mass_width * 2)
    xmax = float(_np.max(xs) + mass_width * 2)
    ax.set_xlim(xmin, xmax)

    yvals = []
    for Larr in (medians, up1, up2, down1, down2):
        mask = ~_np.isnan(Larr)
        if _np.any(mask):
            yvals.extend(Larr[mask].tolist())
    if yvals:
        ymin = min(yvals) * 0.9
        ymax = max(yvals) * 1.1
        if ymin == ymax:
            ymin *= 0.5
            ymax *= 1.5

    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    if show:
        _plt.show()
    return ax


def upper_limit_xsec_br(signal_yield: float,
                        background_yields: Union[List[float], np.ndarray],
                        lumi_ab: float = 1.0,
                        CL: float = 0.95):
    """Median expected 95% CL upper limit on sigma*BR for single-bin counting (uses Poisson CDF + brentq)."""
    from scipy.stats import poisson
    from scipy.optimize import brentq

    signal_yield *= lumi_ab
    background_yields = np.array(background_yields) * lumi_ab
    b = np.sum(background_yields)
    n_obs = int(round(b))

    def cdf_diff(s):
        return poisson.cdf(n_obs, s + b) - (1 - CL)

    s95 = brentq(cdf_diff, 0, 10 * max(1, b + signal_yield))
    xsec_br_limit = s95 / signal_yield if signal_yield > 0 else np.inf
    return s95, xsec_br_limit


def optimize_cut_punzi(hist_signal: Dict[str, Dict[str, Any]],
                       hist_backgrounds: Dict[str, Dict[str, Any]],
                       optimization_params: Dict[str, Any],
                       cut_type: str = "upper",
                       a: float = 2,
                       print_all: bool = False,
                       plot_config: Optional[Dict[str, Any]] = None):
    """Perform cut optimization using Punzi figure-of-merit.

    Returns dict with best_cut, best_estimator, cut_values, estimators.
    """
    histname_to_optimize = optimization_params["histname_to_optimize"]
    cut_values = optimization_params["cut_values"]
    total_signal_ref_hist = optimization_params.get("total_signal_ref_hist", "00_Initial_n_muons")

    estimators = []

    # total initial signal (kept same as notebook: hardcoded factor)
    total_initial_signal = 1 * 0.17
    if total_initial_signal == 0:
        print("Warning: Total initial signal is zero. Efficiencies will be zero.")

    total_sig_no_cut = 0
    for sig_proc, histograms in hist_signal.items():
        hist = histograms.get(histname_to_optimize)
        if hist is None:
            continue
        for i in range(1, hist.GetNbinsX() + 1):
            total_sig_no_cut += hist.GetBinContent(i)
    total_bkg_no_cut = 0
    for bkg_proc, histograms in hist_backgrounds.items():
        hist = histograms.get(histname_to_optimize)
        if hist is None:
            continue
        for i in range(1, hist.GetNbinsX() + 1):
            total_bkg_no_cut += hist.GetBinContent(i)

    base_sig_eff = total_sig_no_cut / total_initial_signal if total_initial_signal > 0 else 0
    base_fom = base_sig_eff / (a / 2 + np.sqrt(total_bkg_no_cut)) if total_bkg_no_cut > 0 else 0
    print(f"Base Punzi FoM with no cut: {base_fom:.6f}")

    for cut in cut_values:
        # background after cut
        total_bkg = 0
        for bkg_proc, histograms in hist_backgrounds.items():
            hist = histograms.get(histname_to_optimize)
            if hist is None:
                continue
            tot_bkg_proc = 0
            for i in range(1, hist.GetNbinsX() + 1):
                bin_center = hist.GetXaxis().GetBinCenter(i)
                if (cut_type == "upper" and bin_center <= cut) or (cut_type == "lower" and bin_center >= cut):
                    tot_bkg_proc += hist.GetBinContent(i)
            total_bkg += tot_bkg_proc

        total_sig = 0
        for sig_proc, histograms in hist_signal.items():
            hist = histograms.get(histname_to_optimize)
            if hist is None:
                continue
            for i in range(1, hist.GetNbinsX() + 1):
                bin_center = hist.GetXaxis().GetBinCenter(i)
                if (cut_type == "upper" and bin_center <= cut) or (cut_type == "lower" and bin_center >= cut):
                    total_sig += hist.GetBinContent(i)

        sig_eff = total_sig / total_initial_signal if total_initial_signal > 0 else 0
        estimated = sig_eff / (a / 2 + np.sqrt(total_bkg)) if total_bkg > 0 else 0
        estimators.append(estimated)
        if print_all:
            print(f"Cut: {cut:.2f}, Signal: {total_sig:.4f}, sig_eff: {sig_eff:.4f}, Bkg: {total_bkg:.4f}, Punzi FoM: {estimated:.6f}")

    if plot_config is None:
        plot_config = {}
    plt.figure(figsize=(10, 6))
    plt.plot(cut_values, estimators, marker="o")
    plt.axhline(y=base_fom, color="r", linestyle="--", label="Base FoM (no cut)")
    plt.legend()
    plt.xlabel(plot_config.get("xlabel", "Cut Value"))
    plt.ylabel(plot_config.get("ylabel", r"Punzi Significance $\epsilon_S / (a/2 + \sqrt{B})$"))
    plt.title(plot_config.get("title", "Cut Optimization"))
    plt.grid(True)
    plt.show()

    best_idx = int(np.argmax(estimators))
    best_cut = cut_values[best_idx]
    best_estimator = estimators[best_idx]
    print(f"Optimal Cut: {best_cut:.3f} -> Max Punzi Significance: {best_estimator:.6f}")

    return {
        "best_cut": best_cut,
        "best_estimator": best_estimator,
        "cut_values": cut_values,
        "estimators": estimators,
    }


def feldman_cousins_upper_limit(n_obs: int, b_expected: float, cl: float = 0.95, mu_step: float = 0.005, mu_max: float = 20.0):
    """Compute upper limit using Feldman-Cousins unified approach (simple scan)."""
    import numpy as _np
    from scipy.stats import poisson

    mus = _np.arange(0, mu_max, mu_step)
    n_outcomes = _np.arange(0, 50)
    last_included_mu = 0.0

    for mu in mus:
        prob_n_given_mu = poisson.pmf(n_outcomes, mu + b_expected)
        mu_best = _np.maximum(0, n_outcomes - b_expected)
        prob_n_given_best = poisson.pmf(n_outcomes, mu_best + b_expected)
        R = _np.divide(prob_n_given_mu, prob_n_given_best, out=_np.zeros_like(prob_n_given_mu), where=prob_n_given_best != 0)
        ranked_indices = _np.argsort(R)[::-1]

        current_sum_prob = 0.0
        is_n_obs_in_region = False
        for idx in ranked_indices:
            n_val = n_outcomes[idx]
            current_sum_prob += prob_n_given_mu[idx]
            if n_val == n_obs:
                is_n_obs_in_region = True
            if current_sum_prob >= cl:
                break

        if is_n_obs_in_region:
            last_included_mu = mu
        else:
            return last_included_mu
    return last_included_mu