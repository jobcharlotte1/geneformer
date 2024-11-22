from kedro.io import AbstractDataset
import os
import pyarrow as pa
import pyarrow.parquet as pq
import json

class GeneformerDataset(AbstractDataset):
    def __init__(self, filepath: str):
        self.filepath = filepath

    def _load(self):
        # Load Arrow data and associated JSON metadata
        arrow_file = os.path.join(self.filepath, "data.arrow")
        metadata_file = os.path.join(self.filepath, "metadata.json")
        table = pq.read_table(arrow_file)
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
        return table, metadata

    def _save(self, data: tuple):
        # Save Arrow table and JSON metadata
        table, metadata = data
        arrow_file = os.path.join(self.filepath, "data.arrow")
        metadata_file = os.path.join(self.filepath, "metadata.json")
        pq.write_table(table, arrow_file)
        with open(metadata_file, "w") as f:
            json.dump(metadata, f)

    def _describe(self):
        return f"Geneformer dataset at {self.filepath}"
