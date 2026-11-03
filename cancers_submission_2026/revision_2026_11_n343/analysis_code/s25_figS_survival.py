"""Final survival figures: full-follow-up fits, five-year KM display and complete sensitivities.
Original plotting implementation remains available in earlier Git history.
"""
from pathlib import Path
import runpy,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"reviewer_followup_20261103/code"))
if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1]/"reviewer_followup_20261103/code/regenerate_survival_figures.py"),run_name="__main__")
