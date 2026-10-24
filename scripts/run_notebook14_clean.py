import json
from pathlib import Path
import nbformat
from nbclient import NotebookClient

p = Path(r"notebooks\14_Surrogate_Model_Communication_Programs.ipynb")
nb = nbformat.read(p, as_version=4)

# The prior cells define X, Y, models and best_models, so a fresh
# kernel cannot safely execute only the final-model and SHAP cells.
# We will instead execute the full notebook to preserve dependencies.
client = NotebookClient(
    nb,
    timeout=900,
    kernel_name="python3",
    allow_errors=False
)

print("Starting full, clean notebook execution...")
client.execute()
nbformat.write(nb, p)
print("Notebook completed and saved successfully.")
