# dataset

Simulated datasets built for a statistical-modelling course, and the analyses of them.

Each round of the exercise works the same way: one person simulates a dataset with a story and a set of
deliberate modelling traps, sends it to a partner with only the README, and the partner analyses it blind.
Each folder below holds one dataset.

| Folder | Dataset | Outcome | Origin |
|---|---|---|---|
| [`project_lumen/`](project_lumen/) | Backyard firefly counts, 25 yards x 20 nights | `fireflies_counted` (count) | Simulated here, sent out for analysis |
| [`project_fulcrum/`](project_fulcrum/) | 500 agonist antibodies against one receptor | `agonism_pct` (continuous, floored at 0) | Received from a classmate, analysed here |

Both datasets are 500 rows with 3-5 predictors and a single outcome.

## project_lumen

Dataset and README as sent to the analysing partner. The generating model and its answer key are
deliberately not published here.

## project_fulcrum

The dataset and README as received, plus the full analysis: a write-up with the recommended model and
held-out performance, one figure, all 500 antibodies ranked by predicted activity, and an inventory of
the traps the dataset contains.

**Note:** `project_fulcrum/analysis/fulcrum_tricky_parts.txt` is effectively an answer key for that
dataset. Read it only after you have analysed the data yourself.
