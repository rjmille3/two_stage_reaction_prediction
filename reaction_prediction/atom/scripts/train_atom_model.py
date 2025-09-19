#!/usr/bin/env python3
# Train a single model with JSON-config hyperparameters and live (per-epoch) plots.
# Class weights are AUTO-BALANCED from training labels.

import os
import json
import time
import gc
import argparse
import h5py
import numpy as np
import sys

# Use a non-interactive backend if no display is available
if os.environ.get("DISPLAY", "") == "":
    import matplotlib
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.regularizers import l2
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint


def parse_args():
    p = argparse.ArgumentParser(
        description="Train a Keras MLP on HDF5 features/targets with live plots. "
                    "Hyperparameters are loaded from a JSON config."
    )
    # Required paths
    p.add_argument("--train_h5", required=True, help="Path to train.hdf5")
    p.add_argument("--val_h5",   required=True, help="Path to val.hdf5")
    p.add_argument("--model_out", required=True, help="Where to save best model (.h5)")
    p.add_argument("--plot_dir",  required=True, help="Directory to save plots & logs")
    p.add_argument("--config",    required=True, help="Path to config.json")

    # Optional: override HDF5 keys if your files differ
    p.add_argument("--features_key", default="features", help="HDF5 dataset key for X")
    p.add_argument("--targets_key",  default="targets",  help="HDF5 dataset key for y")

    return p.parse_args()


def load_hparams(config_path):
    with open(config_path, "r") as f:
        hparams = json.load(f)

    # Basic sanity/defaults in case some fields are missing
    defaults = {
        "num_layers": 6,
        "hidden_dim": 540,
        "dropout": 0.02,
        "l2_reg": 1e-3,
        "lr": 3.28e-05,
        "beta1": 0.866,
        "beta2": 0.998,
        "epochs": 100,
        "batch_size": 512,
        "patience": 100,
        "monitor": "val_loss",   # or "val_accuracy"
        "mode": "min",           # "min" for loss, "max" for accuracy
        "seed": 0
    }
    for k, v in defaults.items():
        hparams.setdefault(k, v)
    return hparams


def load_split(path, features_key="features", targets_key="targets"):
    with h5py.File(path, "r") as f:
        X = f[features_key][:]
        y = np.array(list(f[targets_key][:])).astype(np.float32)
    if y.ndim > 1 and y.shape[1] == 1:
        y = y.reshape(-1)
    return X.astype(np.float32), y


def build_model(input_dim, hparams):
    model = Sequential()
    model.add(Dense(
        hparams["hidden_dim"], activation="relu",
        kernel_regularizer=l2(hparams["l2_reg"]),
        input_shape=(input_dim,)
    ))
    if hparams["dropout"] > 0:
        model.add(Dropout(hparams["dropout"]))

    for _ in range(hparams["num_layers"] - 1):
        model.add(Dense(
            hparams["hidden_dim"], activation="relu",
            kernel_regularizer=l2(hparams["l2_reg"])
        ))
        if hparams["dropout"] > 0:
            model.add(Dropout(hparams["dropout"]))

    model.add(Dense(1, activation="sigmoid"))
    return model


class LivePlotCallback(tf.keras.callbacks.Callback):
    def __init__(self, out_dir, save_prefix="training"):
        super().__init__()
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)
        self.save_prefix = save_prefix
        self.tr_losses, self.va_losses = [], []
        self.tr_accs, self.va_accs = [], []
        self.interactive = (os.environ.get("DISPLAY", "") != "")
        if self.interactive:
            plt.ion()
        # Pre-create figs/axes
        self.fig_loss, self.ax_loss = plt.subplots()
        self.fig_acc, self.ax_acc = plt.subplots()

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        self.tr_losses.append(float(logs.get("loss", np.nan)))
        self.va_losses.append(float(logs.get("val_loss", np.nan)))
        self.tr_accs.append(float(logs.get("accuracy", np.nan)))
        self.va_accs.append(float(logs.get("val_accuracy", np.nan)))
        epochs = np.arange(1, len(self.tr_losses) + 1)

        # Loss plot
        self.ax_loss.clear()
        self.ax_loss.plot(epochs, self.tr_losses, label="train_loss")
        self.ax_loss.plot(epochs, self.va_losses, label="val_loss")
        self.ax_loss.set_title("Loss (Train vs Val)")
        self.ax_loss.set_xlabel("Epoch")
        self.ax_loss.set_ylabel("Binary Crossentropy")
        self.ax_loss.legend()
        self.fig_loss.tight_layout()
        self.fig_loss.canvas.draw()
        if self.interactive:
            self.fig_loss.canvas.flush_events()
        self.fig_loss.savefig(os.path.join(self.out_dir, f"{self.save_prefix}_loss.png"), dpi=150)

        # Accuracy plot
        self.ax_acc.clear()
        self.ax_acc.plot(epochs, self.tr_accs, label="train_acc")
        self.ax_acc.plot(epochs, self.va_accs, label="val_acc")
        self.ax_acc.set_title("Accuracy (Train vs Val)")
        self.ax_acc.set_xlabel("Epoch")
        self.ax_acc.set_ylabel("Accuracy")
        self.ax_acc.legend()
        self.fig_acc.tight_layout()
        self.fig_acc.canvas.draw()
        if self.interactive:
            self.fig_acc.canvas.flush_events()
        self.fig_acc.savefig(os.path.join(self.out_dir, f"{self.save_prefix}_acc.png"), dpi=150)


def main():
    args = parse_args()
    hparams = load_hparams(args.config)

    # Logging setup
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    try:
        tf.get_logger().setLevel('ERROR')
    except Exception:
        pass
    np.random.seed(hparams["seed"])
    tf.random.set_seed(hparams["seed"])

    # Load data
    X_train, y_train = load_split(args.train_h5, args.features_key, args.targets_key)
    X_val,   y_val   = load_split(args.val_h5,   args.features_key, args.targets_key)
    print("Train:", X_train.shape, y_train.shape)
    print("Val:  ", X_val.shape,   y_val.shape, flush=True)

    # --- AUTO-BALANCED class weights from label counts (no fixed weights) ---
    n_pos = int(np.count_nonzero(y_train == 1))
    n_neg = int(np.count_nonzero(y_train == 0))
    total = n_pos + n_neg

    if n_pos == 0 or n_neg == 0:
        raise ValueError(
            f"Dataset imbalance error: one class is missing "
            f"(pos={n_pos}, neg={n_neg}). Cannot compute class weights."
        )
    else:
        # Inverse-frequency normalized so mean weight ≈ 1
        w_neg = (1.0 / n_neg) * (total / 2.0)
        w_pos = (1.0 / n_pos) * (total / 2.0)
        class_weight = {0: w_neg, 1: w_pos}
        print(f"[INFO] Using balanced weights: pos={n_pos}, neg={n_neg} -> {class_weight}")

    # Build model
    K.clear_session()
    gc.collect()

    model = build_model(
        input_dim=X_train.shape[1],
        hparams=hparams
    )

    opt = tf.keras.optimizers.Adam(
        learning_rate=hparams["lr"],
        beta_1=hparams["beta1"],
        beta_2=hparams["beta2"],
        epsilon=1e-8,
        amsgrad=False,
    )

    model.compile(optimizer=opt, loss='binary_crossentropy', metrics=['accuracy'])
    print(model.summary())

    # Callbacks
    os.makedirs(os.path.dirname(args.model_out), exist_ok=True)
    os.makedirs(args.plot_dir, exist_ok=True)

    es = EarlyStopping(
        monitor=hparams["monitor"],
        patience=hparams["patience"],
        mode=hparams["mode"],
        restore_best_weights=True,
        verbose=1
    )

    ckpt = ModelCheckpoint(
        filepath=args.model_out,
        monitor=hparams["monitor"],
        mode=hparams["mode"],
        save_best_only=True,
        verbose=1
    )

    liveplot = LivePlotCallback(out_dir=args.plot_dir, save_prefix="training")

    # Train
    t0 = time.time()
    history = model.fit(
        X_train, y_train,
        validation_data=(X_val, y_val),
        epochs=hparams["epochs"],
        batch_size=hparams["batch_size"],
        shuffle=True,
        class_weight=class_weight,
        callbacks=[es, ckpt, liveplot],
        verbose=1
    )
    print(f"Training time: {time.time()-t0:.2f}s")

    # Save final model (best already saved at model_out)
    base, ext = os.path.splitext(args.model_out)
    final_path = f"{base}_final{ext or '.h5'}"
    model.save(final_path)
    print(f"Saved final model to: {final_path}")

    # Save history and the exact config used
    with open(os.path.join(args.plot_dir, "history.json"), "w") as f:
        json.dump(history.history, f, indent=2)
    with open(os.path.join(args.plot_dir, "config_used.json"), "w") as f:
        json.dump(hparams, f, indent=2)

    # Cleanup
    del model
    K.clear_session()
    gc.collect()


if __name__ == "__main__":
    main()
