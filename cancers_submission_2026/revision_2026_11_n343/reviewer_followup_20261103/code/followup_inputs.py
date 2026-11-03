from pathlib import Path
import sys
FOLLOWUP = Path(__file__).resolve().parents[1]
REVISION = FOLLOWUP.parent
sys.path.insert(0, str(REVISION / "additional_analyses_20261103" / "code"))
from private_inputs import load_survival_observations
def observations():
    a = load_survival_observations()
    import pandas as pd
    public = pd.read_csv(REVISION / "additional_analyses_20261103/inputs/public_annotations_n343.tsv", sep="\t", index_col=0)
    if not a.immune_state.equals(public.immune_state):
        raise ValueError("Restricted input immune_state must equal the canonical 122/144/77 public labels")
    return a
