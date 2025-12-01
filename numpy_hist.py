import numpy as np
from copy import deepcopy


class _Axis:
    def __init__(self, edges: np.ndarray, title: str = ""):
        # edges: array of length nbins+1
        self._edges = np.asarray(edges, dtype=float)
        self._title = str(title) if title is not None else ""

    # ROOT-like API
    def GetTitle(self):
        return self._title

    def SetTitle(self, title: str):
        self._title = str(title)

    def GetNbins(self) -> int:
        return max(len(self._edges) - 1, 0)

    def GetBinLowEdge(self, i: int) -> float:
        nb = self.GetNbins()
        if i < 1:
            i = 1
        if i > nb + 1:
            i = nb + 1
        return float(self._edges[i - 1])

    def GetBinUpEdge(self, i: int) -> float:
        nb = self.GetNbins()
        if i < 1:
            i = 1
        if i > nb:
            i = nb
        return float(self._edges[i])

    def GetBinCenter(self, i: int) -> float:
        low = self.GetBinLowEdge(i)
        up = self.GetBinUpEdge(i)
        return 0.5 * (low + up)

    def FindBin(self, x: float) -> int:
        # Return 1-based bin index; clamp to [1, nbins]
        nb = self.GetNbins()
        if nb <= 0:
            return 1
        # np.searchsorted returns index where to insert to keep order
        idx = int(np.searchsorted(self._edges, x, side="right")) - 1
        # valid interior is [0, nb-1]; translate to [1, nb]
        idx = max(0, min(idx, nb - 1))
        return idx + 1

    # Convenience
    @property
    def edges(self) -> np.ndarray:
        return self._edges


class NumpyHist1D:
    """
    Lightweight numpy-backed 1D histogram wrapper providing a ROOT-TH1 like API used in the notebook.
    Supported methods (1-based indices like ROOT):
      - Clone([new_name]) -> NumpyHist1D
      - Scale(factor) -> self
      - Rebin(factor:int) -> self (in-place)
      - GetTitle(), SetTitle()
      - GetXaxis() -> _Axis, GetYaxis() -> _Axis (dummy for compatibility)
      - GetNbinsX() -> int, GetNbinsY() -> int (always 1)
      - GetEntries() -> float
      - GetBinContent(i:int) -> float, GetBinError(i:int) -> float
    """

    def __init__(self,
                 name: str,
                 title: str,
                 edges: np.ndarray,
                 contents: np.ndarray,
                 errors: np.ndarray | None = None,
                 entries: float | int | None = None,
                 x_title: str = "",
                 y_title: str = ""):
        self._name = str(name) if name is not None else ""
        self._title = str(title) if title is not None else ""
        self._xaxis = _Axis(np.asarray(edges, dtype=float), x_title)
        nbins = self._xaxis.GetNbins()
        c = np.asarray(contents, dtype=float).reshape(-1)
        if c.size != nbins:
            raise ValueError(f"Contents length {c.size} does not match nbins {nbins}.")
        self._contents = c.copy()
        if errors is None:
            # Poisson-like default
            self._errors = np.sqrt(np.clip(self._contents, 0.0, None))
        else:
            e = np.asarray(errors, dtype=float).reshape(-1)
            if e.size != nbins:
                raise ValueError(f"Errors length {e.size} does not match nbins {nbins}.")
            self._errors = e.copy()
        # Entries are number of filled events (unaffected by Scale)
        self._entries = float(entries) if entries is not None else float(np.sum(self._contents))
        # Provide a dummy Y axis for compatibility (used for titles/plots in notebook)
        self._yaxis = _Axis(np.array([0.0, 1.0], dtype=float), y_title)

    # --- Converters ---
    @classmethod
    def from_root_th1(cls, th1) -> "NumpyHist1D":
        # Extract binning, contents, errors, titles from a ROOT.TH1 object
        nb = th1.GetXaxis().GetNbins()
        edges = [th1.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 1)]
        edges.append(th1.GetXaxis().GetBinUpEdge(nb))
        contents = [th1.GetBinContent(i) for i in range(1, nb + 1)]
        errors = [th1.GetBinError(i) for i in range(1, nb + 1)]
        name = th1.GetName() if hasattr(th1, "GetName") else "hist"
        title = th1.GetTitle() if hasattr(th1, "GetTitle") else name
        x_title = th1.GetXaxis().GetTitle() if hasattr(th1.GetXaxis(), "GetTitle") else ""
        y_title = th1.GetYaxis().GetTitle() if hasattr(th1, "GetYaxis") and hasattr(th1.GetYaxis(), "GetTitle") else ""
        entries = th1.GetEntries() if hasattr(th1, "GetEntries") else None
        return cls(name=name,
                   title=title,
                   edges=np.array(edges, dtype=float),
                   contents=np.array(contents, dtype=float),
                   errors=np.array(errors, dtype=float),
                   entries=entries,
                   x_title=x_title,
                   y_title=y_title)

    # --- ROOT-like API ---
    def SetDirectory(self, *_args, **_kwargs):
        # No-op for compatibility
        return None

    def Clone(self, new_name: str | None = None) -> "NumpyHist1D":
        clone = deepcopy(self)
        if new_name is not None:
            clone._name = str(new_name)
        return clone

    def Scale(self, factor: float) -> "NumpyHist1D":
        f = float(factor)
        self._contents *= f
        self._errors *= abs(f)
        # Do not modify entries (mimic ROOT)
        return self

    def Rebin(self, factor: int) -> "NumpyHist1D":
        f = int(factor)
        if f <= 1:
            return self
        nb = self.GetNbinsX()
        if nb % f != 0:
            raise ValueError(f"Rebin factor {f} must divide number of bins {nb}.")
        new_nbins = nb // f
        # New edges by stepping every f bins
        e = self._xaxis.edges
        new_edges = e[::f]
        if new_edges.size != new_nbins + 1:
            # Ensure we include the last edge
            new_edges = np.concatenate([new_edges, e[-1:]])
        # Group-sum contents and errors (errors in quadrature)
        c = self._contents.reshape(new_nbins, f)
        e2 = (self._errors.reshape(new_nbins, f)) ** 2
        new_contents = np.sum(c, axis=1)
        new_errors = np.sqrt(np.sum(e2, axis=1))
        # Replace internals
        self._xaxis = _Axis(new_edges, self._xaxis.GetTitle())
        self._contents = new_contents
        self._errors = new_errors
        # Entries unchanged
        return self

    def GetTitle(self):
        return self._title

    def SetTitle(self, title: str):
        self._title = str(title)

    def GetXaxis(self) -> _Axis:
        return self._xaxis

    def GetYaxis(self) -> _Axis:
        return self._yaxis

    def GetNbinsX(self) -> int:
        return self._xaxis.GetNbins()

    def GetNbinsY(self) -> int:
        return 1

    def GetEntries(self) -> float:
        return self._entries

    def GetBinContent(self, i: int) -> float:
        nb = self.GetNbinsX()
        if i < 1:
            i = 1
        if i > nb:
            i = nb
        return float(self._contents[i - 1])

    def GetBinError(self, i: int) -> float:
        nb = self.GetNbinsX()
        if i < 1:
            i = 1
        if i > nb:
            i = nb
        return float(self._errors[i - 1])
    
    # Get 
    def Merge(self, other: "NumpyHist1D") -> "NumpyHist1D":
        if not np.array_equal(self._xaxis.edges, other._xaxis.edges):
            raise ValueError("Cannot merge histograms with different binning.")
        new_contents = self._contents + other._contents
        new_errors = np.sqrt(self._errors**2 + other._errors**2)
        new_entries = self._entries + other._entries
        merged_hist = NumpyHist1D(
            name=self._name,
            title=self._title,
            edges=self._xaxis.edges,
            contents=new_contents,
            errors=new_errors,
            entries=new_entries,
            x_title=self._xaxis.GetTitle(),
            y_title=self._yaxis.GetTitle()
        )
        return merged_hist

    # Convenience accessors
    @property
    def name(self) -> str:
        return self._name

    @property
    def title(self) -> str:
        return self._title

    @property
    def edges(self) -> np.ndarray:
        return self._xaxis.edges

    @property
    def contents(self) -> np.ndarray:
        return self._contents

    @property
    def errors(self) -> np.ndarray:
        return self._errors


class NumpyHist2D:
    """
    Numpy-backed 2D histogram with a ROOT-TH2-like API used in the notebook.
    Supported methods:
      - Clone([new_name]) -> NumpyHist2D
      - Scale(factor) -> self
      - Rebin(factor:int | tuple[int,int]) -> self (in-place)
      - GetTitle(), SetTitle()
      - GetXaxis() -> _Axis, GetYaxis() -> _Axis
      - GetNbinsX() -> int, GetNbinsY() -> int
      - GetEntries() -> float
      - GetBinContent(ix:int, iy:int) -> float, GetBinError(ix:int, iy:int) -> float
    Indices are 1-based like ROOT.
    """

    def __init__(self,
                 name: str,
                 title: str,
                 x_edges: np.ndarray,
                 y_edges: np.ndarray,
                 contents: np.ndarray,
                 errors: np.ndarray | None = None,
                 entries: float | int | None = None,
                 x_title: str = "",
                 y_title: str = "",
                 z_title: str = ""):
        self._name = str(name) if name is not None else ""
        self._title = str(title) if title is not None else ""
        self._xaxis = _Axis(np.asarray(x_edges, dtype=float), x_title)
        self._yaxis = _Axis(np.asarray(y_edges, dtype=float), y_title)
        nbx = self._xaxis.GetNbins()
        nby = self._yaxis.GetNbins()
        c = np.asarray(contents, dtype=float)
        if c.shape != (nbx, nby):
            raise ValueError(f"Contents shape {c.shape} does not match (nbx, nby)=({nbx}, {nby}).")
        self._contents = c.copy()
        if errors is None:
            self._errors = np.sqrt(np.clip(self._contents, 0.0, None))
        else:
            e = np.asarray(errors, dtype=float)
            if e.shape != (nbx, nby):
                raise ValueError(f"Errors shape {e.shape} does not match (nbx, nby)=({nbx}, {nby}).")
            self._errors = e.copy()
        self._entries = float(entries) if entries is not None else float(np.sum(self._contents))
        # Keep a dummy Z title in y-axis object if needed; ROOT would use TH1::GetYaxis for vertical axis and colorbar title via hist title; we store separately if needed
        self._z_title = str(z_title) if z_title is not None else ""

    # --- Converters ---
    @classmethod
    def from_root_th2(cls, th2) -> "NumpyHist2D":
        nbx = th2.GetXaxis().GetNbins()
        nby = th2.GetYaxis().GetNbins()
        x_edges = [th2.GetXaxis().GetBinLowEdge(i) for i in range(1, nbx + 1)]
        x_edges.append(th2.GetXaxis().GetBinUpEdge(nbx))
        y_edges = [th2.GetYaxis().GetBinLowEdge(j) for j in range(1, nby + 1)]
        y_edges.append(th2.GetYaxis().GetBinUpEdge(nby))
        contents = np.zeros((nbx, nby), dtype=float)
        errors = np.zeros((nbx, nby), dtype=float)
        for i in range(1, nbx + 1):
            for j in range(1, nby + 1):
                contents[i - 1, j - 1] = th2.GetBinContent(i, j)
                errors[i - 1, j - 1] = th2.GetBinError(i, j)
        name = th2.GetName() if hasattr(th2, "GetName") else "hist2d"
        title = th2.GetTitle() if hasattr(th2, "GetTitle") else name
        x_title = th2.GetXaxis().GetTitle() if hasattr(th2.GetXaxis(), "GetTitle") else ""
        y_title = th2.GetYaxis().GetTitle() if hasattr(th2.GetYaxis(), "GetTitle") else ""
        # No distinct Z axis in TH2; use title or empty
        entries = th2.GetEntries() if hasattr(th2, "GetEntries") else None
        return cls(name=name,
                   title=title,
                   x_edges=np.array(x_edges, dtype=float),
                   y_edges=np.array(y_edges, dtype=float),
                   contents=contents,
                   errors=errors,
                   entries=entries,
                   x_title=x_title,
                   y_title=y_title,
                   z_title="")

    # --- ROOT-like API ---
    def SetDirectory(self, *_args, **_kwargs):
        return None

    def Clone(self, new_name: str | None = None) -> "NumpyHist2D":
        clone = deepcopy(self)
        if new_name is not None:
            clone._name = str(new_name)
        return clone

    def Scale(self, factor: float) -> "NumpyHist2D":
        f = float(factor)
        self._contents *= f
        self._errors *= abs(f)
        return self

    def Rebin(self, factor) -> "NumpyHist2D":
        # factor can be int or (fx, fy)
        if isinstance(factor, (list, tuple, np.ndarray)):
            if len(factor) != 2:
                raise ValueError("For 2D Rebin, provide a single int or a tuple (fx, fy).")
            fx, fy = int(factor[0]), int(factor[1])
        else:
            fx = fy = int(factor)
        if fx <= 1 and fy <= 1:
            return self
        nbx = self.GetNbinsX()
        nby = self.GetNbinsY()
        if fx <= 0 or fy <= 0:
            raise ValueError("Rebin factors must be positive integers.")
        if nbx % fx != 0 or nby % fy != 0:
            raise ValueError(f"Rebin factors (fx={fx}, fy={fy}) must divide (nbx={nbx}, nby={nby}).")
        new_nbx = nbx // fx
        new_nby = nby // fy

        # New edges
        ex = self._xaxis.edges
        ey = self._yaxis.edges
        new_ex = ex[::fx]
        if new_ex.size != new_nbx + 1:
            new_ex = np.concatenate([new_ex, ex[-1:]])
        new_ey = ey[::fy]
        if new_ey.size != new_nby + 1:
            new_ey = np.concatenate([new_ey, ey[-1:]])

        # Group-sum contents and errors (quadrature) by block reshaping
        c = self._contents.reshape(new_nbx, fx, new_nby, fy).swapaxes(1, 2)
        e2 = (self._errors.reshape(new_nbx, fx, new_nby, fy).swapaxes(1, 2)) ** 2
        new_contents = c.sum(axis=(2, 3))  # sum over fx, fy
        new_errors = np.sqrt(e2.sum(axis=(2, 3)))

        # Replace
        self._xaxis = _Axis(new_ex, self._xaxis.GetTitle())
        self._yaxis = _Axis(new_ey, self._yaxis.GetTitle())
        self._contents = new_contents
        self._errors = new_errors
        return self

    def GetTitle(self):
        return self._title

    def SetTitle(self, title: str):
        self._title = str(title)

    def GetXaxis(self) -> _Axis:
        return self._xaxis

    def GetYaxis(self) -> _Axis:
        return self._yaxis

    def GetNbinsX(self) -> int:
        return self._xaxis.GetNbins()

    def GetNbinsY(self) -> int:
        return self._yaxis.GetNbins()

    def GetEntries(self) -> float:
        return self._entries

    def GetBinContent(self, ix: int, iy: int) -> float:
        nbx = self.GetNbinsX()
        nby = self.GetNbinsY()
        if ix < 1:
            ix = 1
        if iy < 1:
            iy = 1
        if ix > nbx:
            ix = nbx
        if iy > nby:
            iy = nby
        return float(self._contents[ix - 1, iy - 1])

    def GetBinError(self, ix: int, iy: int) -> float:
        nbx = self.GetNbinsX()
        nby = self.GetNbinsY()
        if ix < 1:
            ix = 1
        if iy < 1:
            iy = 1
        if ix > nbx:
            ix = nbx
        if iy > nby:
            iy = nby
        return float(self._errors[ix - 1, iy - 1])
    
    def Merge(self, other: "NumpyHist2D") -> "NumpyHist2D":
        if not np.array_equal(self._xaxis.edges, other._xaxis.edges):
            raise ValueError("Cannot merge histograms with different X binning.")
        if not np.array_equal(self._yaxis.edges, other._yaxis.edges):
            raise ValueError("Cannot merge histograms with different Y binning.")
        new_contents = self._contents + other._contents
        new_errors = np.sqrt(self._errors**2 + other._errors**2)
        new_entries = self._entries + other._entries
        merged_hist = NumpyHist2D(
            name=self._name,
            title=self._title,
            x_edges=self._xaxis.edges,
            y_edges=self._yaxis.edges,
            contents=new_contents,
            errors=new_errors,
            entries=new_entries,
            x_title=self._xaxis.GetTitle(),
            y_title=self._yaxis.GetTitle(),
            z_title=self._z_title
        )
        return merged_hist

    # Convenience
    @property
    def name(self) -> str:
        return self._name

    @property
    def title(self) -> str:
        return self._title

    @property
    def x_edges(self) -> np.ndarray:
        return self._xaxis.edges

    @property
    def y_edges(self) -> np.ndarray:
        return self._yaxis.edges

    @property
    def contents(self) -> np.ndarray:
        return self._contents

    @property
    def errors(self) -> np.ndarray:
        return self._errors

