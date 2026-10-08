"""Run LimiX-2 zero-shot classification in an isolated process.

Why a separate process: TabFM (JAX) and LimiX-2 (PyTorch) each initialize their own
threading/OpenMP runtimes. After a large JAX inference session, torch inference in
the SAME process slows down by orders of magnitude (measured on this machine:
37 s standalone vs >40 min in-process). Process isolation keeps both honest.

Usage: python run_limix_asd.py <in.npz> <out.npz> [config_name]
  in.npz:  X_train (n_tr, p), y_train (n_tr,), X_test (n_te, p)  [float arrays]
  out.npz: proba (n_te, n_classes), t_load (s), t_pred (s)
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# LimiX's distributed utilities expect these even in single-process runs
for _k, _v in [("RANK", "0"), ("WORLD_SIZE", "1"),
               ("MASTER_ADDR", "127.0.0.1"), ("MASTER_PORT", "29500")]:
    os.environ.setdefault(_k, _v)

import numpy as np
import torch
from huggingface_hub import hf_hub_download

LIMIX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "LimiX")
sys.path.insert(0, LIMIX_DIR)
from inference.predictor import LimiXPredictor  # noqa: E402


def main():
    in_path, out_path = sys.argv[1], sys.argv[2]
    config_name = sys.argv[3] if len(sys.argv) > 3 else "cls_course_trimmed6_v2.json"

    data = np.load(in_path)
    ckpt = hf_hub_download(repo_id="stableai-org/LimiX-2", filename="LimiX-2.ckpt")

    t0 = time.time()
    clf = LimiXPredictor(
        device=torch.device("cpu"),
        model_path=ckpt,
        inference_config=os.path.join(LIMIX_DIR, "config", config_name),
    )
    t_load = time.time() - t0

    t0 = time.time()
    proba = clf.predict(data["X_train"], data["y_train"], data["X_test"],
                        task_type="Classification")
    t_pred = time.time() - t0

    np.savez(out_path, proba=proba, t_load=t_load, t_pred=t_pred)
    print(f"limix done: load {t_load:.1f}s, predict {t_pred:.1f}s")


if __name__ == "__main__":
    main()
