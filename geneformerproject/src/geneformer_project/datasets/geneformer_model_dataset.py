import os
from pathlib import Path
import shutil
import json
import pickle
from kedro.io.core import AbstractDataSet, DataSetError

class GeneformerModelDataset(AbstractDataSet):
    def __init__(self, filepath: str):
        """Dataset for handling Geneformer classifier save format.
        
        Args:
            filepath (str): Path to the directory where the model is saved or loaded.
        """
        self._filepath = Path(filepath)
        self._required_files = [
            "checkpoints/",
            "predictions.pkl",
            "config.json",
            "model.safetensors",
            "training_args.bin",
        ]

    def _validate_files(self, path: Path):
        """Ensure the required files exist in the directory."""
        for file in self._required_files:
            if not (path / file).exists():
                raise DataSetError(f"Missing required file or folder: {file}")

    def _load(self) -> dict:
        """Load the Geneformer model save directory as a dictionary."""
        self._validate_files(self._filepath)
        result = {
            "checkpoints": list((self._filepath / "checkpoints").glob("**/*")),
            "predictions": pickle.load(open(self._filepath / "predictions.pkl", "rb")),
            "config": json.load(open(self._filepath / "config.json", "r")),
            "model": self._filepath / "model.safetensors",
            "training_args": self._filepath / "training_args.bin",
        }
        return result

    def _save(self, data: dict) -> None:
        """Save the Geneformer model directory."""
        self._filepath.mkdir(parents=True, exist_ok=True)
        
        # Save checkpoints (assumes list of paths)
        checkpoints_dir = self._filepath / "checkpoints"
        checkpoints_dir.mkdir(exist_ok=True)
        for checkpoint_file in data["checkpoints"]:
            shutil.copy(checkpoint_file, checkpoints_dir)
        
        # Save predictions
        with open(self._filepath / "predictions.pkl", "wb") as f:
            pickle.dump(data["predictions"], f)
        
        # Save config
        with open(self._filepath / "config.json", "w") as f:
            json.dump(data["config"], f)
        
        # Save model and training arguments
        shutil.copy(data["model"], self._filepath / "model.safetensors")
        shutil.copy(data["training_args"], self._filepath / "training_args.bin")

    def _describe(self):
        return dict(filepath=str(self._filepath))
