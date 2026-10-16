"""Result 1, step 03: spacing test for several sets at once, with and without AP-1 sites in repeat sequence
(code in lib/spacing_followup.py).

Usage: python 03_spacing_sets_and_repeat_masking.py --set NAME RRE.bed CONTROL.bed GENOME.2bit GENOME.sizes
       [--set ...] --out result.json
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from spacing_followup import main  # noqa: E402

if __name__ == "__main__":
    main()
