"""Result 1, step 02: primary AP-1 spacing test for one enhancer set (code in lib/spacing_core.py).

Usage: python 02_ap1_spacing_test.py --rre RRE.bed --ctrl CONTROL_POOL.bed --genome-2bit GENOME.2bit
       --genome-sizes GENOME.sizes --out spacing.json
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from spacing_core import main  # noqa: E402

if __name__ == "__main__":
    main()
