#!/usr/bin/env python3
"""Run InstaNovo de novo prediction with a torch>=2.6 weights_only workaround.

InstaNovo (<=1.2.2) hardcodes torch.load(weights_only=True) at
transformer/model.py:153. torch>=2.6 changed the weights_only default to True and
its restricted unpickler rejects the v1.2.0 checkpoint's collections.defaultdict
(the TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD env var is ignored when the argument is passed
explicitly). The InstaDeepAI checkpoints come from a trusted source, so we allowlist
defaultdict and force weights_only=False, then dispatch to the InstaNovo CLI unchanged.

Arguments pass through verbatim, e.g.:
    instanovo_predict.py transformer predict -d spectra.mzML -o preds.csv --denovo
"""

import collections
import sys

import torch

torch.serialization.add_safe_globals([collections.defaultdict])
_orig_load = torch.load
torch.load = lambda *a, **k: _orig_load(*a, **{**k, "weights_only": False})

from instanovo.cli import instanovo_entrypoint  # noqa: E402

if __name__ == "__main__":
    sys.exit(instanovo_entrypoint())
