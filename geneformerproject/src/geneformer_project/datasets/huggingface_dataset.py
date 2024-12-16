from kedro.io import AbstractDataset

class HuggingFaceModelDataset(AbstractDataset):
    def __init__(self, filepath: str):
        self._filepath = filepath

    def _load(self):
        from transformers import AutoModel
        return AutoModel.from_pretrained(self._filepath)

    def _save(self, model):
        model.save_pretrained(self._filepath)

    def _describe(self):
        return dict(filepath=self._filepath)