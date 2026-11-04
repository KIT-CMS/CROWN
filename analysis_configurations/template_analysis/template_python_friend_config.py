"""
Template for a python (ML) friend config, evaluated by KingMaker's PythonFriend task.

Unlike CROWN friend configs, this one is not compiled: the jobs import it and call
build_config(), load_model() and evaluate() in the ML environment (by default an LCG
view, see ml_setup_script in KingMaker_luigi.cfg). Any model type works, as long as
evaluate() returns one numpy array per output.

The jobs run on GPUs with htcondor_request_gpus = 1 in [PythonFriend], so write
load_model()/evaluate() to use a GPU if there is one and the CPU otherwise:
onnxruntime as below, tensorflow does this on its own, torch e.g. with
    device = "cuda" if torch.cuda.is_available() else "cpu"
and model.to(device), torch.as_tensor(..., device=device), result.cpu().numpy().

Keep module level imports light: the config is also imported on the submission host
to check the inputs, where torch/onnxruntime/... are not available.

The models used here can be created with e.g.
    import torch
    model = torch.nn.Sequential(torch.nn.Linear(6, 16), torch.nn.ReLU(), torch.nn.Linear(16, 1), torch.nn.Sigmoid())
    torch.onnx.export(model, torch.zeros(1, 6), "payloads/python_friend_example/nominal.onnx",
                      input_names=["x"], output_names=["y"], dynamic_axes={"x": {0: "n"}, "y": {0: "n"}})
"""

from __future__ import annotations

import os

# names the model was trained with, in training order
MODEL_INPUTS = ["mu1_pt", "mu2_pt", "mu1_eta", "mu2_eta", "pair_mass", "pair_pt"]


def build_config(era, sample_type, scope, shifts):
    """
    Returns which models are evaluated on which columns.

    evaluations: one entry per model. "shift" is "nominal" or the name of a shift,
        outputs of shift models are written as <output>__<shift>. Only shifts selected
        via --shifts are evaluated (same matching as CROWN), nominal always.
    input_map: model input name -> ntuple column, so the model does not have to be
        retrained if the ntuples name a quantity differently. Can also be given per
        evaluation, e.g. to point a shift model to shifted columns.
    empty_values: model input -> constant used if the column is missing. Only for
        tests on incomplete ntuples, needs --allow-empty-values in production.
    input_ranges, input_groups: only for the sensitivity test (--sensitivity-test),
        see the README. input_ranges needs an entry for every model input.
    """
    return {
        "evaluations": [
            {
                "shift": "nominal",
                "model": "payloads/python_friend_example/nominal.onnx",
                "inputs": MODEL_INPUTS,
                "outputs": ["example_score"],
            },
            {
                "shift": "MuonIDUp",
                "model": "payloads/python_friend_example/MuonIDUp.onnx",
                "inputs": MODEL_INPUTS,
                "outputs": ["example_score"],
                # e.g. if this model was trained on shifted inputs:
                # "input_map": {"mu1_pt": "pt_1__MuonIDUp"},
            },
        ],
        "input_map": {
            "mu1_pt": "pt_1",
            "mu2_pt": "pt_2",
            "mu1_eta": "eta_1",
            "mu2_eta": "eta_2",
            "pair_mass": "mm_pair_mass",
            "pair_pt": "mm_pair_pt",
        },
        # "empty_values": {"pair_pt": -10.0},
        # random inputs of the sensitivity test, keyed by model input name
        "input_ranges": {
            "mu1_pt": {"min": 20.0, "max": 200.0},  # uniform
            "mu2_pt": {"min": 10.0, "max": 100.0},
            "mu1_eta": {"mean": 0.0, "std": 1.2},  # gaussian
            "mu2_eta": {"values": [-2.0, -1.0, 0.0, 1.0, 2.0], "weights": [1, 2, 4, 2, 1]},
            "pair_mass": {"min": 70.0, "max": 110.0},
            # with an empty_values entry for pair_pt, the empty value can be mixed in:
            # "pair_pt": {"min": 0.0, "max": 200.0, "empty_fraction": 0.1},
            "pair_pt": pair_pt,  # derived from inputs sampled before
        },
        "input_groups": {"muon_pts": ["mu1_pt", "mu2_pt"]},
        "batch_size": 100_000,
    }


def pair_pt(rng, n, data):
    """
    Derived input of the sensitivity test: random opening angle between the two muons.
    """
    import numpy as np

    dphi = rng.uniform(-np.pi, np.pi, n)
    return np.sqrt(
        data["mu1_pt"] ** 2 + data["mu2_pt"] ** 2 + 2 * data["mu1_pt"] * data["mu2_pt"] * np.cos(dphi)
    )


def load_model(path, evaluation):
    """
    Called once per job and evaluation with the local path of the model file.
    """
    import onnxruntime

    options = onnxruntime.SessionOptions()
    # stay within the cores requested from the batch system
    options.intra_op_num_threads = int(os.environ.get("PYTHON_FRIEND_NTHREADS", 1))
    options.inter_op_num_threads = 1
    # the GPU if the view has CUDA and the job got one, the CPU otherwise
    providers = [
        provider
        for provider in ("CUDAExecutionProvider", "CPUExecutionProvider")
        if provider in onnxruntime.get_available_providers()
    ]
    return onnxruntime.InferenceSession(path, options, providers=providers)


def evaluate(model, data, evaluation):
    """
    data: model input name -> numpy array of one batch of events.
    Returns output name -> numpy array with one value per event.
    """
    import numpy as np

    x = np.stack([data[name] for name in evaluation["inputs"]], axis=1)
    y = model.run(None, {model.get_inputs()[0].name: x.astype(np.float32)})[0]
    return {"example_score": y[:, 0]}
