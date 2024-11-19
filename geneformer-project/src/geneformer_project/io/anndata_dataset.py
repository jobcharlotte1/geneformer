from kedro.io import AbstractDataSet
from anndata import AnnData, read_h5ad

class AnnDataSet(AbstractDataSet):
    """
    Custom Kedro dataset to handle AnnData (.h5ad) files.
    """

    def __init__(self, filepath: str):
        """
        Initializes the AnnDataSet.
        Args:
            filepath (str): Path to the .h5ad file.
        """
        self._filepath = filepath

    def _load(self) -> AnnData:
        """
        Load the AnnData object from the .h5ad file.
        Returns:
            AnnData: The loaded AnnData object.
        """
        return read_h5ad(self._filepath)

    def _save(self, data: AnnData) -> None:
        """
        Save an AnnData object to an .h5ad file.
        Args:
            data (AnnData): The AnnData object to save.
        """
        data.write(self._filepath)

    def _describe(self) -> dict:
        """
        Describe the dataset.
        Returns:
            dict: A dictionary describing the dataset.
        """
        return dict(filepath=self._filepath)
