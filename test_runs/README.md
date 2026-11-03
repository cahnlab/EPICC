# test_runs/

Scratch space for local pipeline test/integration runs. **Gitignored** —
nothing here is tracked (this README is the only exception, force-added).

Point a test run's outputs here so the working tree stays clean:

```bash
epicc run \
  --profile <your-site-profile> \
  --options tests/integration/data/test_options_colcen.yaml \
  --samples tests/integration/data/test_samples_colcen.tsv \
  --output-dir test_runs/colcen/results \
  --genome-dir test_runs/colcen/genomes \
  2>&1 | tee test_runs/colcen/run.log
```

Delete a run with `rm -rf test_runs/<name>`. Older runs used top-level
`results_test_*/` and `genomes_test_*/` dirs (also gitignored); migrate
those into here when convenient.
