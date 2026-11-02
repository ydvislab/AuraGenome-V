"""Write the synthetic example used in the README and the tests (the paper's data cannot be redistributed).

A region of 600 substitutions and 120 tracks in 12 families (biosamples), each family sharing a latent effect pattern.
The cell type of interest ("our cells") is not profiled; its true effects follow family 4 with a different scale and
offset, plus measurement noise. Family 4's biosample is "fibroblast of lung", so --anchor lung names the right
family (--anchor fibroblast would also take the dermal fibroblasts); --anchor keratinocyte names a related but wrong
one (family 5, correlated 0.5 with family 4).

usage: python examples/make_example.py   ->  examples/example_predictions.tsv.gz, examples/example_measurements.tsv
"""
import os

import numpy as np
import pandas as pd

rng = np.random.default_rng(1)
n, n_fam, per = 600, 12, 10
names = ["fibroblast of dermis", "hepatocyte", "K562", "B cell", "fibroblast of lung", "keratinocyte", "astrocyte",
         "HUVEC", "monocyte", "neural progenitor", "myotube", "osteoblast"]
heavy = rng.standard_t(3, size=(n, n_fam))                       # heavy-tailed, like variant effects
lat = heavy.copy()
lat[:, 5] = 0.5 * lat[:, 4] + np.sqrt(1 - 0.25) * lat[:, 5]       # keratinocyte partly shares the lung-fibroblast signal
ids = [f"chr7:{1_000_000 + i // 3}:{'ACGT'[i % 4]}>{'CGTA'[i % 4]}" for i in range(n)]
rows = []
for f in range(n_fam):
    for r in range(per):
        assay = ["DNASE", "CAGE", "ATAC"][r % 3]
        z = lat[:, f] * rng.uniform(0.5, 2) + rng.normal(0, 0.3, n)
        rows.append(pd.DataFrame(dict(variant_id=ids, output_type=assay, track_name=f"{assay} {names[f]} rep{r}",
                                      track_strand=".", biosample_name=names[f], raw_score=0.05 * z)))
pred = pd.concat(rows, ignore_index=True)
truth = 0.3 + 0.8 * lat[:, 4] + rng.normal(0, 0.5, n)
here = os.path.dirname(os.path.abspath(__file__))
pred.to_csv(os.path.join(here, "example_predictions.tsv.gz"), sep="\t", index=False, float_format="%.5g")
pd.DataFrame(dict(variant_id=ids, effect=truth)).to_csv(os.path.join(here, "example_measurements.tsv"), sep="\t",
                                                        index=False, float_format="%.5g")
print(f"wrote {len(pred)} prediction rows and {n} measurements to {here}")
