"""Central, validated path registry for the whole project.

Single source of truth for every directory. Import `P` and address artifacts
by intent, never by hand-built strings:

    from paths import P
    df = pd.read_parquet(P.interim / "GSE16879_expr.parquet")
    fig.savefig(P.out("02_clock") / "cv.png")          # machine output, disposable
    # promotion (final artifact) is a deliberate act -> see P.promote_figure()

Two-tier artifact model (do not blur these):
  outputs/  = whatever a script emits automatically. Regenerable, disposable,
              machine-owned, never hand-edited. Keyed by script stem via P.out().
  results/  = the curated subset you consciously PROMOTE as paper-final.
              figures/ + tables/. A copy lives here only after promotion.

Layout is resolved relative to this file (repo_root = parent of src/), so the
tree is portable: no absolute paths baked in.
"""
from __future__ import annotations

import shutil
from pathlib import Path

from pydantic import BaseModel, ConfigDict, computed_field


class Paths(BaseModel):
    """Validated directory registry. Frozen: paths are constants, not state."""

    model_config = ConfigDict(frozen=True)

    root: Path = Path(__file__).resolve().parent.parent

    # --- data tiers ---
    @computed_field
    @property
    def data(self) -> Path:
        return self.root / "data"

    @computed_field
    @property
    def raw(self) -> Path:
        """Immutable, read-only source data (GWAS, eQTLGen, GEO, GTEx)."""
        return self.data / "raw"

    @computed_field
    @property
    def external(self) -> Path:
        """Reference resources (gene sets, annotations)."""
        return self.data / "external"

    @computed_field
    @property
    def interim(self) -> Path:
        """Derived, regenerable data (scored tables, parquet matrices)."""
        return self.data / "interim"

    # --- code & narrative ---
    @computed_field
    @property
    def src(self) -> Path:
        return self.root / "src"

    @computed_field
    @property
    def notebooks(self) -> Path:
        return self.root / "notebooks"

    @computed_field
    @property
    def journal(self) -> Path:
        return self.root / "research_journal"

    @computed_field
    @property
    def paper(self) -> Path:
        return self.root / "paper"

    # --- artifact tiers ---
    @computed_field
    @property
    def outputs(self) -> Path:
        """Machine auto-emit, disposable. Use P.out(stem) for a script subdir."""
        return self.root / "outputs"

    @computed_field
    @property
    def results(self) -> Path:
        """Curated, paper-final. Populated only via promotion."""
        return self.root / "results"

    @computed_field
    @property
    def figures(self) -> Path:
        return self.results / "figures"

    @computed_field
    @property
    def tables(self) -> Path:
        return self.results / "tables"

    # --- helpers ---
    def out(self, script_stem: str) -> Path:
        """Per-script output dir (created on demand), e.g. P.out("02_clock")."""
        d = self.outputs / script_stem
        d.mkdir(parents=True, exist_ok=True)
        return d

    def promote_figure(self, src_file: Path) -> Path:
        """Copy a machine output into results/figures/ — the 'this is final' act."""
        self.figures.mkdir(parents=True, exist_ok=True)
        dst = self.figures / Path(src_file).name
        shutil.copy2(src_file, dst)
        return dst

    def promote_table(self, src_file: Path) -> Path:
        """Copy a machine output into results/tables/ — the 'this is final' act."""
        self.tables.mkdir(parents=True, exist_ok=True)
        dst = self.tables / Path(src_file).name
        shutil.copy2(src_file, dst)
        return dst


P = Paths()

if __name__ == "__main__":
    for k, v in P.model_dump().items():
        print(f"{k:10} {v}")
