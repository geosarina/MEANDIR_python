"""Reconstruct Table S2 from the supplement PDF using character coordinates,
so the scenario columns (which interleave in naive text extraction) are
assigned correctly. Prints a clean grid: 18 (variable, end-member) rows x 9
columns (All/Mainstem/Slough x scenario 1/2/3).
"""
import sys
import types

# stub cryptography (its rust binding panics here; only needed for encrypted PDFs)
for n in ["cryptography", "cryptography.hazmat", "cryptography.hazmat.primitives",
          "cryptography.hazmat.primitives.ciphers", "cryptography.hazmat.backends"]:
    sys.modules[n] = types.ModuleType(n)
ci = sys.modules["cryptography.hazmat.primitives.ciphers"]
ci.Cipher = ci.algorithms = ci.modes = object
sys.modules["cryptography.hazmat.backends"].default_backend = lambda *a, **k: None

from pdfminer.high_level import extract_pages
from pdfminer.layout import LTTextContainer, LTTextLine, LTChar, LTTextLineHorizontal

PDF = "/root/.claude/uploads/feec57dc-e16b-557c-afd4-47f9242d36db/80bf6e88-Kemeney_Supplement_2022gb007644sup0001supporting_information_sis01.pdf"
PAGE = 10  # 0-based index of the Table S2 page

ROWS = [  # (variable, end-member) in table order
    ("DIC", "Carbonate"), ("DIC", "Corg oxidation"), ("DIC", "Degassing"),
    ("Ca", "Carbonate"), ("Ca", "Evaporite"), ("Ca", "Silicate"), ("Ca", "Precipitation"),
    ("Mg", "Carbonate"), ("Mg", "Silicate"), ("Mg", "Precipitation"),
    ("Na", "Silicate"), ("Na", "Precipitation"),
    ("K", "Silicate"), ("K", "Precipitation"),
    ("Cl", "Precipitation"),
    ("SO4", "H2SO4 production"), ("SO4", "Evaporite"), ("SO4", "Precipitation"),
]
COLS = ["All|1", "All|2", "All|3", "Main|1", "Main|2", "Main|3", "Slough|1", "Slough|2", "Slough|3"]


def _chars(page):
    """recursively yield LTChar objects."""
    stack = list(page)
    while stack:
        el = stack.pop()
        if isinstance(el, LTChar):
            yield el
        elif isinstance(el, LTTextContainer) or hasattr(el, "__iter__"):
            try:
                stack.extend(list(el))
            except TypeError:
                pass


def tokens():
    """Build value-tokens from individual characters, splitting on x-gaps so
    that two numbers sharing a baseline become two tokens."""
    out = []
    for page in extract_pages(PDF, page_numbers=[PAGE]):
        chars = [c for c in _chars(page) if c.get_text().strip()]
        # group by row (y), then split into tokens on x gaps
        chars.sort(key=lambda c: (-round((c.y0 + c.y1) / 2 / 3), c.x0))
        rows = {}
        for c in chars:
            yk = round((c.y0 + c.y1) / 2 / 3)
            rows.setdefault(yk, []).append(c)
        for yk, cs in rows.items():
            cs.sort(key=lambda c: c.x0)
            cur, x0, x1 = "", None, None
            for c in cs:
                if x1 is not None and c.x0 - x1 > 2.0:  # gap -> new token
                    _emit(out, cur, x0, x1, yk * 3)
                    cur, x0 = "", None
                if x0 is None:
                    x0 = c.x0
                cur += c.get_text()
                x1 = c.x1
            _emit(out, cur, x0, x1, yk * 3)
    return out


def _emit(out, text, x0, x1, y):
    tt = text.replace("%", "").replace("−", "-").strip()
    if tt == "-" or _isnum(tt):
        out.append((tt, (x0 + x1) / 2, y))


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def cluster(vals, tol):
    vals = sorted(vals)
    groups = [[vals[0]]]
    for v in vals[1:]:
        if v - groups[-1][-1] <= tol:
            groups[-1].append(v)
        else:
            groups.append([v])
    return [sum(g) / len(g) for g in groups]


def main():
    toks = tokens()
    print(f"value tokens found on page {PAGE}: {len(toks)}")
    ys = cluster([t[2] for t in toks], tol=4)
    xs = cluster([t[1] for t in toks], tol=12)
    ys = sorted(ys, reverse=True)  # top to bottom
    xs = sorted(xs)                # left to right
    print(f"row clusters: {len(ys)}  col clusters: {len(xs)}")

    # build grid
    grid = {}
    for txt, x, y in toks:
        ri = min(range(len(ys)), key=lambda i: abs(ys[i] - y))
        ci_ = min(range(len(xs)), key=lambda i: abs(xs[i] - x))
        grid[(ri, ci_)] = txt
    # dump raw grid for inspection
    print("\nraw grid (row index : values left->right):")
    for ri in range(len(ys)):
        cells = [grid.get((ri, ci_), ".") for ci_ in range(len(xs))]
        print(f"{ri:2d}: " + "  ".join(f"{c:>6s}" for c in cells))


if __name__ == "__main__":
    main()
